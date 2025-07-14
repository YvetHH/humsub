from io import StringIO, BytesIO
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
import torch
from torch import nn
import esm
from esm.pretrained import load_model_and_alphabet_hub
from Bio import SeqIO

from snakemake.shell import shell
import logging, traceback
import os, glob
from collections import defaultdict

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

def read_sequence_from_fasta(file_path, sequence_id):
    for record in SeqIO.parse(file_path, "fasta"):
        if record.id == sequence_id:
            return str(record.seq)
    return None

def calculate_median_seq_len(file_path, sequences_id):
    seq_len_list = []
    for seq_id in sequences_id:
        seq = read_sequence_from_fasta(file_path, seq_id)
        len_seq = len(seq)
        seq_len_list.append(len_seq)
    return np.median(seq_len_list)


models_path = snakemake.config['model_path']

model_names = [
    "ESMForSingleMutationPosOuter",
    "ESMForSingleMutationPosConcat",
    "ESMForSingleMutation_pos_cat_cls",
    "ESMForSingleMutation_pos",
    "ESMForSingleMutation_cls",
]

model_names = [f"{models_path}/{model_name}" for model_name in model_names]

import torch.nn.functional as F

HIDDEN_UNITS_POS_CONTACT = 5


class ESMForSingleMutationPosConcat(nn.Module):
    def __init__(self):
        super().__init__()
        self.esm2, _ = esm.pretrained.esm2_t33_650M_UR50D()
        self.fc1 = nn.Linear(1280 * 2, HIDDEN_UNITS_POS_CONTACT)
        self.fc2 = nn.Linear(HIDDEN_UNITS_POS_CONTACT, 1)

    def forward(self, token_ids1, token_ids2, pos):
        outputs1 = self.esm2.forward(token_ids1, repr_layers=[33])["representations"][33]
        outputs2 = self.esm2.forward(token_ids2, repr_layers=[33])["representations"][33]
        outputs1_pos = outputs1[:, pos + 1]
        outputs2_pos = outputs2[:, pos + 1]
        outputs_pos_concat = torch.cat((outputs1_pos, outputs2_pos), 2)
        fc1_outputs = F.relu(self.fc1(outputs_pos_concat))
        logits = self.fc2(fc1_outputs)
        return logits


HIDDEN_UNITS_POS_OUTER = 5


class ESMForSingleMutationPosOuter(nn.Module):
    def __init__(self):
        super().__init__()
        self.esm2, _ = esm.pretrained.esm2_t33_650M_UR50D()
        self._freeze_esm2_layers()
        self.fc1 = nn.Linear(1280 * 1280, HIDDEN_UNITS_POS_OUTER)
        self.fc2 = nn.Linear(HIDDEN_UNITS_POS_OUTER, 1)

    def _freeze_esm2_layers(self):
        total_blocks = 33
        initial_layers = 2
        layers_per_block = 16
        num_freeze_blocks = total_blocks - 3
        for _, param in list(self.esm2.named_parameters())[: initial_layers + layers_per_block * num_freeze_blocks]:
            param.requires_grad = False

    def forward(self, token_ids1, token_ids2, pos):
        outputs1 = self.esm2.forward(token_ids1, repr_layers=[33])["representations"][33]
        outputs2 = self.esm2.forward(token_ids2, repr_layers=[33])["representations"][33]
        outputs1_pos = outputs1[:, pos + 1]
        outputs2_pos = outputs2[:, pos + 1]
        outer_prod = outputs1_pos.unsqueeze(3) @ outputs2_pos.unsqueeze(2)
        outer_prod_view = outer_prod.view(outer_prod.shape[0], outer_prod.shape[1], -1)
        fc1_outputs = F.relu(self.fc1(outer_prod_view))
        logits = self.fc2(fc1_outputs)
        return logits


class ESMForSingleMutation_pos(nn.Module):
    def __init__(self):
        super().__init__()
        self.esm1v, self.esm1v_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.classifier = nn.Linear(1280, 1)
        self.const1 = torch.nn.Parameter(torch.ones((1, 1280)))
        self.const2 = torch.nn.Parameter(-1 * torch.ones((1, 1280)))

    def forward(self, token_ids1, token_ids2, pos):
        outputs1 = self.esm1v.forward(token_ids1, repr_layers=[33])["representations"][33]
        outputs2 = self.esm1v.forward(token_ids2, repr_layers=[33])["representations"][33]
        outputs = self.const1 * outputs1[:, pos + 1, :] + self.const2 * outputs2[:, pos + 1, :]
        logits = self.classifier(outputs)
        return logits


class ESMForSingleMutation_cls(nn.Module):
    def __init__(self):
        super().__init__()
        self.esm1v, self.esm1v_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.classifier = nn.Linear(1280, 1)
        self.const1 = torch.nn.Parameter(torch.ones((1, 1280)))
        self.const2 = torch.nn.Parameter(-1 * torch.ones((1, 1280)))

    def forward(self, token_ids1, token_ids2, pos):
        outputs1 = self.esm1v.forward(token_ids1, repr_layers=[33])["representations"][33]
        outputs2 = self.esm1v.forward(token_ids2, repr_layers=[33])["representations"][33]
        outputs = self.const1 * outputs1[:, 0, :] + self.const2 * outputs2[:, 0, :]
        logits = self.classifier(outputs.unsqueeze(0))
        return logits


class ESMForSingleMutation_pos_cat_cls(nn.Module):
    def __init__(self):
        super().__init__()
        self.esm1v, self.esm1v_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.classifier = nn.Linear(1280 * 2, 1)
        self.const1 = torch.nn.Parameter(torch.ones((1, 1280)))
        self.const2 = torch.nn.Parameter(-1 * torch.ones((1, 1280)))

    def forward(self, token_ids1, token_ids2, pos):
        outputs1 = self.esm1v.forward(token_ids1, repr_layers=[33])["representations"][33]
        outputs2 = self.esm1v.forward(token_ids2, repr_layers=[33])["representations"][33]
        cls_out = self.const1 * outputs1[:, 0, :] + self.const2 * outputs2[:, 0, :]
        pos_out = self.const1 * outputs1[:, pos + 1, :] + self.const2 * outputs2[:, pos + 1, :]
        outputs = torch.cat([cls_out.unsqueeze(0), pos_out], axis=-1)
        logits = self.classifier(outputs)
        return logits

def analyze_function(function, df):
    logging.info(f"Analyzing function {function} . . .")
    output_df = pd.DataFrame()
    function_df = df.query('function == @function')
    effect_dict = defaultdict()

    for pos in function_df['Position'].unique():
        tmp_function_df = function_df.query("Position == @pos")
        reference_group = tmp_function_df['Group'].iloc[0]
        reference_gene = tmp_function_df['Reference Sequence ID'].iloc[0]
        reference_alel = tmp_function_df.query("Group == @reference_group")['Variant'].iloc[0]
        other_aleles = tmp_function_df.query("Group != @reference_group")['Variant'].tolist()

        for other_alel in set(other_aleles):
            if reference_alel == other_alel:
                continue
            mutation_code = f"{reference_alel}{pos}{other_alel}"
            seq = read_sequence_from_fasta(f"output/genes/{snakemake.wildcards.cluster95}.faa", reference_gene)

            # Get wildtype sequence, mutation position and mutated sequence
            wt_aa = mutation_code[0]
            mut_aa = mutation_code[-1]
            mut_pos = int(mutation_code[1:-1]) - 1

            wt = seq
            tt = list(seq)
            tt[mut_pos] = mut_aa
            mut = "".join(tt)

            model = torch.load(f"{models_path}/ESMForSingleMutation_cls", map_location=torch.device("cpu"))
            esm2_alphabet = model.esm1v_alphabet
            esm2batch_converter = esm2_alphabet.get_batch_converter()
            _, _, esm2_batch_tokens1 = esm2batch_converter([("", wt[:1022])])
            _, _, esm2_batch_tokens2 = esm2batch_converter([("", mut[:1022])])
            esm2_batch_tokens1 = esm2_batch_tokens1.cpu()
            esm2_batch_tokens2 = esm2_batch_tokens2.cpu()

            res = []
            for model_name in model_names:
                model = torch.load(model_name, map_location=torch.device("cpu"))
                model.eval()
                model.cpu() # HERE

                with torch.no_grad():
                    res.append(
                        model(token_ids1=esm2_batch_tokens1, token_ids2=esm2_batch_tokens2, pos=torch.LongTensor([mut_pos]))
                        .cpu()
                        .numpy()
                    )
            res = np.mean(res)
            output_df = pd.concat([output_df, pd.DataFrame({"position":pos, 'reference_group':reference_group,
                                                            "reference_gene":reference_gene, "reference_allel":reference_alel, 'other_alel':other_alel,'effect':res}, index=[function])], axis=0)
    logging.info(f"The function {function} is successfully analyzed!")
    return output_df  


def run_analyze_function(function, df):
    try:
        return analyze_function(function, df)
    except Exception as e:
        logging.info(f"Error processing function {function}: {str(e)}")
        return pd.DataFrame()  # Return an empty DataFrame in case of an error



msa_batch = pd.read_csv(snakemake.input.batch_df, sep='\t').drop("Unnamed: 0",axis=1)
msa_batch['Group'] = msa_batch['Group'].astype(str).str.rjust(10,'0')

existing_subspecies = msa_batch['Group'].unique()
os.makedirs(f"output/prostata/tmp/batch_output/{snakemake.wildcards.cluster95}", exist_ok=True)

final_df = pd.DataFrame()
for function in msa_batch['function'].unique():
    function_df = run_analyze_function(function, msa_batch)
    final_df = pd.concat([final_df,function_df],axis=0)


final_df.to_csv(snakemake.output.variant_effect, sep='\t')
