---
Title: Two-Stage vs Single-Stage Training Analysis
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
    - Path: SATURN/data/gene_embeddings.py
      Note: |-
        Protein embedding integration - loading mechanism for precomputed embeddings
        Protein embedding integration mechanism - loading precomputed embeddings
    - Path: SATURN/losses/triplet_margin_loss.py
      Note: Metric learning loss formulation - triplet margin loss computation
    - Path: SATURN/miners/triplet_margin_miner.py
      Note: |-
        Triplet mining implementation - cross-species constraint enforcement
        Triplet mining implementation details - cross-species constraint enforcement
    - Path: SATURN/model/saturn_model.py
      Note: |-
        Core model implementations analyzed - SATURNPretrainModel architecture (encoder-decoder)
        Core model architectures analyzed - SATURNPretrainModel (encoder-decoder)
    - Path: SATURN/train-saturn.py
      Note: |-
        Training pipeline implementation - pretrain_saturn() function details pretraining loop with ZINB+L1+protein similarity losses
        Training pipeline implementation details - pretrain_saturn() function shows pretraining loop with ZINB+L1+protein similarity losses
    - Path: SATURN/utils/loss_and_miner_utils.py
      Note: |-
        Triplet mining utilities - get_species_triplet_indices() cross-species implementation
        Triplet mining utility functions - get_species_triplet_indices() cross-species implementation
ExternalSources: []
Summary: Technical analysis of two-stage (pretraining + metric learning) versus single-stage end-to-end training for SATURN architecture
LastUpdated: 2025-12-04T08:32:59.451455546-05:00
---



# Two-Stage vs Single-Stage Training Analysis

## Goal

This document provides a technical analysis comparing two-stage training (pretraining followed by metric learning) versus single-stage end-to-end training approaches for the SATURN architecture. The analysis evaluates architectural requirements, objective conflicts, parameter convergence, and training stability.

## Context

SATURN employs a two-stage training pipeline: (1) pretraining with a conditional autoencoder to learn macrogene representations, and (2) metric learning with triplet loss for cross-species alignment. This analysis examines whether this separation is architecturally necessary or if a unified end-to-end approach could achieve similar or superior results.

## Technical Analysis

### Current Two-Stage Architecture

**Stage 1: Pretraining**
- Model: `SATURNPretrainModel` (conditional autoencoder)
- Input: Gene expression counts (log-transformed) + species encoding + optional batch labels
- Architecture: Genes → Macrogenes → Embeddings → Reconstructed Genes
- Loss: `L_pretrain = L_ZINB + λ_l1 * ||p_weights||_1 + λ_pe * L_pe_sim`
  - `L_ZINB`: Zero-Inflated Negative Binomial reconstruction loss
  - `L1`: Sparsity regularization on macrogene weights
  - `L_pe_sim`: Protein embedding similarity loss
- Duration: ~200 epochs (configurable)
- Output: Frozen macrogene weights, macrogene values per cell

**Stage 2: Metric Learning**
- Model: `SATURNMetricModel` (encoder only, weights copied from pretrain)
- Input: Macrogene values (frozen from pretraining)
- Architecture: Macrogenes → Embeddings (no decoder)
- Loss: `L_metric = max(0, margin + d(anchor, positive) - d(anchor, negative))`
- Triplet mining: Cross-species pairs, optional MNN, species balancing
- Output: Normalized embeddings for cross-species alignment

### Key Technical Considerations

#### 1. Objective Conflict Resolution

**Problem**: Reconstruction and metric learning objectives operate on different optimization landscapes:
- Reconstruction loss (`L_ZINB`) optimizes macrogenes to capture gene expression patterns
- Metric learning loss optimizes embeddings to align cells across species
- These objectives are complementary but require sequential optimization

**Evidence**: The codebase shows macrogene weights are frozen after pretraining (unless `--unfreeze_macrogenes`). The `SATURNMetricModel` discards the decoder entirely, indicating reconstruction is not needed during metric learning.

**Analysis**: Simultaneous optimization risks gradient conflicts. The reconstruction objective pulls macrogenes toward expression pattern fidelity, while metric learning pulls toward cross-species alignment. Two-stage training allows macrogenes to converge to a stable representation before alignment begins.

#### 2. Macrogene Stabilization

**Initialization**: Macrogenes are initialized via K-means clustering on protein embeddings (default: 2000 centroids). Score functions convert distances to weights: `default` (`log(1 + 1/rank)^2 * 2`), `one_hot`, or `smoothed`.

**Stabilization Requirement**: Macrogenes must learn meaningful functional structure before metric learning. The reconstruction task forces macrogenes to capture biologically meaningful patterns. Without pretraining, macrogene weights would be pulled in conflicting directions—learning functional roles while simultaneously optimizing alignment.

**Evidence**: `make_centroids()` in `model/saturn_model.py` performs initialization. `SATURNPretrainModel` initializes `p_weights` from centroid scores and learns them during pretraining. Macrogene values are extracted post-pretraining for metric learning.

**Analysis**: Pretraining provides a dedicated phase for macrogene structure learning. Single-stage training risks macrogene collapse or failure to converge to interpretable functional modules.

#### 3. ZINB Parameter Convergence

**Complexity**: ZINB loss requires learning multiple species-specific parameters:
- Mean: `μ = exp(library_size) * softmax(macrogene_scale)`
- Dispersion: `θ = exp(px_r)` (species-specific `px_rs` dictionary)
- Zero-inflation: `π = sigmoid(px_dropout)` (species-specific decoders)

**Convergence Requirements**: Library size normalization (`library = log(sum(gene_counts))`), species-specific dispersion parameters, and zero-inflation probabilities must converge for accurate count modeling.

**Evidence**: `get_reconstruction_loss()` uses `ZeroInflatedNegativeBinomial` from scvi-tools. `px_rate`, `px_r`, `px_drop` computed in `SATURNPretrainModel.forward()`. Species-specific `px_dropout_decoders` and `px_rs` are ParameterDict structures.

**Analysis**: ZINB parameter learning requires focused optimization. Single-cell data exhibits species-specific overdispersion and zero-inflation patterns. Combining ZINB reconstruction with metric learning risks poor parameter convergence, propagating errors into alignment quality.

#### 4. Regularization and Collapse Prevention

**L1 Regularization**: `l1_penalty * ||p_weights||_1` encourages sparsity, preventing macrogenes from becoming redundant.

**Protein Embedding Similarity**: `pe_sim_penalty * MSE(cosine_sim(macrogene_embeds), cosine_sim(protein_embeds))` maintains alignment with functional protein space.

**Collapse Risk**: Without pretraining, metric learning alone could cause macrogenes to collapse into a few dominant patterns, losing functional diversity. The L1 and protein similarity losses prevent this during pretraining.

**Analysis**: Pretraining provides a regularization phase that prevents macrogene collapse. Single-stage training would require careful loss weighting to balance alignment with regularization, increasing hyperparameter sensitivity.

#### 5. Transferability and Modularity

**Transferability**: Pretrained macrogene weights can be frozen and reused across experiments. Macrogene values are extracted and used as input to `SATURNMetricModel`, enabling modular reuse.

**Modularity**: The separation allows pretrained macrogenes to be transferred without retraining. The encoder-decoder architecture is only needed during pretraining; metric learning uses encoder-only.

**Analysis**: Two-stage training enables model reuse and transfer learning. Single-stage training loses this modularity, requiring full retraining for each experiment.

#### 6. Protein Embedding Refinement

**Initialization Signal**: Protein embeddings provide functional similarity signal independent of gene names. K-means clustering on stacked embeddings initializes macrogenes in functionally meaningful space.

**Refinement Process**: Pretraining refines protein-guided initialization. Protein embeddings provide the seed; pretraining cultivates it into a macrogene representation adapted to expression data while maintaining protein-based functional relationships.

**Evidence**: `data/gene_embeddings.py` loads precomputed embeddings (ESM1b default: 1280-dim). `p_weights_embeddings` layer in `SATURNPretrainModel` encodes macrogene similarity. Protein embedding similarity loss maintains alignment during pretraining.

**Analysis**: Single-stage training would pull macrogenes between protein similarity (initialization), reconstruction (ZINB), and alignment (metric learning) simultaneously. This multi-objective conflict could prevent coherent functional module learning.

### Single-Stage Alternative: Feasibility Analysis

**Potential Benefits**:
- End-to-end optimization could allow macrogenes to adapt directly to metric learning objective
- Potential efficiency gains from unified optimization
- Simplified training pipeline

**Implementation Challenges**:
- Loss weighting: Requires careful balancing of reconstruction and metric losses
- Curriculum learning: Would need gradual transition from reconstruction to alignment (effectively soft two-stage)
- Parameter convergence: ZINB parameters may not converge with competing objectives
- Hyperparameter sensitivity: Increased complexity in loss balancing

**Assessment**: A single-stage approach with curriculum learning (gradual loss weighting) is conceptually similar to two-stage training but less stable. The current two-stage approach provides explicit separation of concerns and proven convergence properties.

### Recommendations

**Maintain Two-Stage Training**:
1. **Architectural Necessity**: Objective conflicts require sequential optimization
2. **Stability**: Pretraining prevents macrogene collapse and ensures parameter convergence
3. **Interpretability**: Two-stage preserves macrogene functional meaning
4. **Transferability**: Enables model reuse and modularity

**Potential Enhancement**: Consider `--unfreeze_macrogenes` during metric learning to allow fine-tuning of macrogene weights while maintaining pretrained initialization. This provides a middle ground between fully frozen and fully trainable macrogenes.

**Future Exploration**: Single-stage training with adaptive loss weighting and curriculum learning deserves experimental validation, but current evidence strongly supports two-stage as the safer, more stable approach.

## Implementation Details

### Pretraining Loss Components

```python
# ZINB Reconstruction Loss
L_ZINB = -log P(x | μ, θ, π)
where:
  μ = exp(library_size) * softmax(macrogene_scale)
  θ = exp(px_r)  # species-specific
  π = sigmoid(px_dropout)  # species-specific

# L1 Regularization
L_l1 = l1_penalty * ||p_weights||_1

# Protein Embedding Similarity
L_pe_sim = pe_sim_penalty * MSE(
    cosine_sim(macrogene_embeds),
    cosine_sim(protein_embeds)
)

# Total Pretraining Loss
L_pretrain = L_ZINB + L_l1 + L_pe_sim
```

### Metric Learning Loss

```python
# Triplet Margin Loss
L_metric = max(0, margin + d(anchor, positive) - d(anchor, negative))

# Distance metric (default: cosine)
d(x, y) = 1 - (x · y) / (||x|| * ||y||)

# Constraints:
# - anchor and positive from different species
# - optional MNN for positive pair selection
# - optional species balancing
```

### Training Pipeline

```python
# Stage 1: Pretraining
pretrain_model = SATURNPretrainModel(...)
pretrain_model = pretrain_saturn(
    pretrain_model, pretrain_loader, optimizer,
    device, nepochs=200, ...
)

# Extract macrogene values
macrogene_values = extract_macrogenes(pretrain_model, data)

# Stage 2: Metric Learning
metric_model = SATURNMetricModel(...)
metric_model.encoder.load_state_dict(pretrain_model.encoder.state_dict())
train(metric_model, triplet_loss, miner, ...)
```

## Related

- [SATURN Architecture and Implementation Analysis](./01-saturn-architecture-and-implementation-analysis.md) - Comprehensive architecture reference
- `SATURN/model/saturn_model.py` - Model implementations
- `SATURN/train-saturn.py` - Training pipeline
- `SATURN/losses/triplet_margin_loss.py` - Metric learning loss
