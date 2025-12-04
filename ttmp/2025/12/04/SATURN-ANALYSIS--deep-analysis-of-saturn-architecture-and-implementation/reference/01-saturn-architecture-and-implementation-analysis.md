---
Title: SATURN Architecture and Implementation Analysis
Ticket: SATURN-ANALYSIS
Status: active
Topics:
    - machine-learning
    - embeddings
    - metric-learning
DocType: reference
Intent: long-term
Owners: []
RelatedFiles:
    - Path: ../../../../../../SATURN/data/gene_embeddings.py
      Note: Protein embedding loading and gene filtering logic
    - Path: ../../../../../../SATURN/data/multi_species_data.py
      Note: Multi-species dataset classes and collation
    - Path: ../../../../../../SATURN/losses/base_metric_loss_function.py
      Note: Base class for metric learning losses
    - Path: ../../../../../../SATURN/miners/triplet_margin_miner.py
      Note: Triplet mining strategy for cross-species alignment
    - Path: ../../../../../../SATURN/model/saturn_model.py
      Note: Core model implementations (SATURNPretrainModel
    - Path: ../../../../../../SATURN/pretrain_utils.py
      Note: Pretraining utilities and KLD scheduling
    - Path: ../../../../../../SATURN/score_adata.py
      Note: Evaluation metrics and scoring functions
    - Path: ../../../../../../SATURN/train-saturn.py
      Note: Main training pipeline orchestrating all stages
ExternalSources: []
Summary: 'Comprehensive technical analysis of SATURN: architecture, control flow, mathematical foundations, package structure, and implementation details'
LastUpdated: 2025-12-04T08:26:29.629652721-05:00
---



# SATURN Architecture and Implementation Analysis

## Goal

This document provides a comprehensive technical analysis of SATURN (Single-cell Analysis Tool for Unifying RNA-seq across species using Neural networks), covering its architecture, control flow, mathematical foundations, package structure, and implementation details. It serves as a reference for understanding how SATURN integrates single-cell RNA-seq datasets across species using protein embeddings and metric learning.

## Context

SATURN is a deep learning approach that couples gene expression data with protein representations learned from large protein language models (PLMs) for cross-species integration. The key innovation is mapping cells from all datasets to a shared space of functionally related genes called "macrogenes." SATURN uniquely detects functionally related genes co-expressed across species by leveraging protein sequence similarity.

**Key Components:**
- Macrogene initialization using K-means clustering on protein embeddings
- Pretraining with a conditional autoencoder (scVI-style ZINB loss)
- Fine-tuning cell clusters with weakly supervised metric learning

## Architecture Overview

### High-Level Pipeline

SATURN operates in three main stages:

1. **Macrogene Initialization**: K-means clustering on protein embeddings to create initial macrogene weights
2. **Pretraining**: Conditional autoencoder learns to reconstruct gene expression from macrogenes
3. **Metric Learning**: Triplet-based metric learning aligns cells across species in embedding space

### Core Models

#### SATURNPretrainModel

The pretraining model (`model/saturn_model.py`) implements a conditional autoencoder:

**Architecture:**
- **Input**: Gene expression counts (log-transformed) + species one-hot encoding + optional batch labels
- **Macrogene Layer**: Linear transformation from genes to macrogenes using learnable weights `p_weights`
- **Encoder**: Multi-layer MLP (hidden_dim → embed_dim) with optional VAE reparameterization
- **Decoder**: Conditional decoder that takes embedding + species/batch info → macrogene scale → gene expression

**Key Components:**
```python
# Gene to macrogene mapping (learnable)
self.p_weights = nn.Parameter(gene_scores.float().t().log())

# Encoder: macrogenes → embeddings
self.encoder = nn.Sequential(
    full_block(self.num_cl, hidden_dim, dropout),
    full_block(hidden_dim, embed_dim, dropout)
)

# Decoder: embeddings + species → macrogenes → genes
self.px_decoder = nn.Sequential(
    full_block(embed_dim + num_species + num_batch_labels, hidden_dim, dropout)
)
```

**Loss Function:**
- **Reconstruction Loss**: Zero-Inflated Negative Binomial (ZINB) loss for count data
- **L1 Regularization**: `l1_penalty * ||p_weights||_1` to encourage sparsity
- **Protein Embedding Similarity Loss**: `pe_sim_penalty * MSE(cosine_sim(macrogene_embeds), cosine_sim(protein_embeds))` to align macrogene embeddings with protein embedding space

#### SATURNMetricModel

The metric learning model operates on macrogene values (frozen from pretraining):

**Architecture:**
- **Input**: Macrogene values (from pretrained model)
- **Encoder**: Same encoder architecture as pretrain model (weights copied)
- **Output**: Normalized embeddings for metric learning

**Key Difference**: No decoder - only encodes macrogenes to embeddings for distance computation.

### Control Flow

#### Training Pipeline (`train-saturn.py`)

**Phase 1: Data Preparation**
1. Load AnnData objects for each species
2. Filter genes to those with protein embeddings
3. Select highly variable genes (default: 8000) using Seurat v3 method
4. Create species-specific gene index mappings

**Phase 2: Macrogene Initialization**
1. Stack protein embeddings from all species
2. Run K-means clustering (default: 2000 centroids)
3. Convert distances to scores using ranking function:
   - `default`: `log(1 + 1/rank)^2 * 2`
   - `one_hot`: Binary assignment to closest centroid
   - `smoothed`: `1/rank`
4. Initialize `p_weights` with log-transformed scores

**Phase 3: Pretraining**
1. Create `SATURNPretrainModel` with initialized macrogene weights
2. Train with conditional autoencoder loss:
   - Forward pass: genes → macrogenes → embeddings → reconstructed genes
   - Loss: ZINB reconstruction + L1 + protein embedding similarity
   - Optional: Balance loss by cell type abundance
3. Freeze macrogene weights (unless `--unfreeze_macrogenes`)

**Phase 4: Metric Learning**
1. Extract macrogene values for all cells using pretrained model
2. Create `SATURNMetricModel` (encoder copied from pretrain)
3. Train with triplet margin loss:
   - Mine triplets: (anchor, positive, negative) where positive/negative are from different species
   - Loss: `max(0, margin + d(anchor, positive) - d(anchor, negative))`
   - Optional: Use mutual nearest neighbors (MNN) for mining
   - Optional: Equalize triplet sampling across species

**Phase 5: Evaluation**
1. Extract final embeddings
2. Score using cross-species classification metrics
3. Output AnnData with embeddings and macrogene values

## Mathematical Foundations

### Macrogene Concept

Macrogenes represent functionally related gene clusters. The mapping from genes to macrogenes is:

```
macrogene_j = Σ_i (w_ij * gene_i)
```

where `w_ij` is the weight from gene `i` to macrogene `j`, initialized from protein embedding similarity.

### Pretraining Loss

**Reconstruction Loss (ZINB):**
```
L_recon = -log P(x | μ, θ, π)
```

where:
- `μ = exp(library_size) * softmax(macrogene_scale)` (mean)
- `θ = exp(px_r)` (dispersion)
- `π = sigmoid(px_dropout)` (zero-inflation probability)

**Regularization:**
```
L_total = L_recon + λ_l1 * ||W||_1 + λ_pe * L_pe_sim
```

where `L_pe_sim` encourages macrogene embeddings to preserve protein embedding similarity structure.

### Metric Learning Loss

**Triplet Margin Loss:**
```
L_triplet = max(0, margin + d(a, p) - d(a, n))
```

where:
- `a` = anchor embedding
- `p` = positive embedding (same cell type, different species)
- `n` = negative embedding (different cell type)
- `d` = cosine distance

**Triplet Mining:**
- **Cross-species**: Ensures anchor and positive are from different species
- **Semihard**: Negative is further than positive but still violates margin
- **MNN**: Uses mutual nearest neighbors to find better positive pairs

### Distance Metrics

**Cosine Similarity** (default):
```
d(x, y) = 1 - (x · y) / (||x|| * ||y||)
```

**Other options**: Lp distance, dot product similarity, SNR distance

## Package Structure

### Core Modules

```
SATURN/
├── model/
│   └── saturn_model.py          # SATURNPretrainModel, SATURNMetricModel
├── data/
│   ├── gene_embeddings.py        # Load protein embeddings
│   └── multi_species_data.py     # Multi-species dataset classes
├── losses/                       # Metric learning losses
│   ├── base_metric_loss_function.py
│   ├── triplet_margin_loss.py
│   └── ... (25+ loss functions)
├── miners/                       # Triplet/pair mining strategies
│   ├── base_miner.py
│   ├── triplet_margin_miner.py
│   └── ... (12+ miners)
├── distances/                    # Distance metrics
│   ├── cosine_similarity.py
│   ├── lp_distance.py
│   └── ...
├── reducers/                     # Loss reduction strategies
│   ├── mean_reducer.py
│   ├── threshold_reducer.py
│   └── ...
├── testers/                      # Evaluation utilities
│   ├── base_tester.py
│   ├── global_embedding_space.py
│   └── ...
├── utils/                        # Helper functions
│   ├── accuracy_calculator.py
│   ├── inference.py
│   └── ...
├── train-saturn.py              # Main training script
├── pretrain_utils.py             # Pretraining utilities
└── score_adata.py                # Evaluation/scoring functions
```

### Key Data Structures

**ExperimentDatasetMulti**: Multi-species dataset that:
- Stores gene expression per species
- Handles species-specific indexing
- Supports batch labels for batch correction

**ExperimentDatasetMultiEqual**: Equalizes cell counts across species by:
- Sampling with replacement from smaller species
- Ensures balanced training across species

**Collate Function**: `multi_species_collate_fn` groups batches by species, returning:
```python
{
    'species1': (data_tensor, labels_tensor, ref_labels_tensor, batch_labels_tensor),
    'species2': (...),
    ...
}
```

## Implementation Details

### Protein Embedding Integration

**Loading (`data/gene_embeddings.py`):**
1. Maps gene symbols to protein IDs
2. Loads precomputed protein embeddings (ESM1b, MSA1b, protXL, ESM2)
3. Filters genes to those with embeddings in all species
4. Returns stacked embedding tensor: `(num_genes_total, embedding_dim)`

**Supported Models:**
- ESM1b (default): 1280-dim embeddings
- ESM2: Updated ESM model
- MSA1b: Multiple sequence alignment model
- protXL: ProtTrans-XL model

### Macrogene Initialization

**K-means Clustering:**
```python
kmeans = KMeans(n_clusters=num_macrogenes, random_state=seed)
kmeans.fit(protein_embeddings)
distances = kmeans.transform(embeddings)  # (num_genes, num_macrogenes)
```

**Score Conversion:**
- Rank distances per gene (rank 1 = closest centroid)
- Convert ranks to scores using chosen function
- Store as initial `p_weights` (log-transformed)

### Pretraining Details

**Conditional Decoding:**
- Species one-hot encoding concatenated with embeddings
- Optional batch labels for batch correction
- Species-specific decoders for gene expression reconstruction

**VAE Option:**
- Optional variational bottleneck: `z ~ N(μ, σ²)` where `μ, σ = encoder(macrogenes)`
- Reparameterization trick: `z = μ + σ * ε`, `ε ~ N(0,1)`
- KL divergence regularization (cycled weight)

**Balancing:**
- Optional cell type balancing: weight loss by inverse cell type frequency
- Prevents rare cell types from being ignored

### Metric Learning Details

**Triplet Mining (`miners/triplet_margin_miner.py`):**
- **Cross-species mining**: Ensures anchor and positive are from different species
- **MNN option**: Uses mutual nearest neighbors to find better positive pairs
- **Ref labels**: Optional coarse labels for additional supervision

**Triplet Types:**
- **Hard**: Negative closer to anchor than positive
- **Semihard**: Negative further than positive but violates margin
- **Easy**: Already satisfies margin (filtered out)

**Species Balancing:**
- Optional equalization of triplet counts per species
- Prevents bias toward larger species

### Evaluation Metrics (`score_adata.py`)

**Cross-Species Classification:**
1. Train logistic regression on species 1 embeddings
2. Test on species 2 embeddings
3. Map predictions using cell type mapping
4. Compute accuracy metrics

**Additional Metrics:**
- **Reannotation**: Louvain clustering + majority voting
- **KNN Alignment**: Cross-species nearest neighbor accuracy
- **Centroid Matching**: Cluster centroid nearest neighbor matching

## Key Design Decisions

### Why Macrogenes?

1. **Functional Alignment**: Protein embeddings capture functional similarity better than gene names
2. **Cross-Species Compatibility**: Orthologs may have different names but similar sequences
3. **Dimensionality Reduction**: Reduces from ~8000 genes to ~2000 macrogenes
4. **Interpretability**: Macrogenes represent functional gene modules

### Why Two-Stage Training?

1. **Stability**: Pretraining learns macrogene structure before metric learning
2. **Regularization**: Pretraining loss prevents macrogenes from collapsing
3. **Transferability**: Pretrained macrogenes can be reused across experiments

### Why Triplet Loss?

1. **Weak Supervision**: Only needs cell type labels, not explicit alignments
2. **Cross-Species Focus**: Naturally handles different species
3. **Flexibility**: Can incorporate MNN and other mining strategies

## Dependencies

**Core:**
- PyTorch (1.10.2+)
- scvi-tools (ZINB loss)
- scanpy (single-cell analysis)
- anndata (data structures)

**ML Libraries:**
- scikit-learn (K-means, classification)
- pytorch-metric-learning (losses, miners, distances)

**Data:**
- pandas, numpy
- scipy

## Usage Examples

### Basic Training

```python
python train-saturn.py \
    --in_data data.csv \
    --in_label_col cell_type \
    --ref_label_col CL_class_coarse \
    --num_macrogenes 2000 \
    --hv_genes 8000 \
    --pretrain_epochs 200 \
    --epochs 50 \
    --embedding_model ESM1b
```

### With Batch Correction

```python
python train-saturn.py \
    --in_data data.csv \
    --non_species_batch_col tissue \
    --use_batch_labels \
    ...
```

### Scoring Results

```python
python score_adata.py \
    --adata results.h5ad \
    --species1 human \
    --species2 mouse \
    --ct_map_path cell_type_mapping.csv \
    --label labels2
```

## Related

- SATURN paper: https://www.biorxiv.org/content/10.1101/2023.02.03.526939v1
- scvi-tools: https://scvi-tools.org/
- pytorch-metric-learning: https://kevinmusgrave.github.io/pytorch-metric-learning/
