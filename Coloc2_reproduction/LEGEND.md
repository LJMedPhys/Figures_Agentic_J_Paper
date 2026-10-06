**Fig. X | Agentic-J reproduces a published Coloc 2 colocalization result.**
Wide-field images of *Drosophila* S2 cells stained for Su(Hw), phosphorylated H2Av (γH2Av) and DAPI, in isotonic medium
(6 fields) or after 250 mM NaCl osmotic stress (11 fields), from Simmons et al. 2022 (*PLoS Genet*,
doi:10.1371/journal.pgen.1010396, Fig. 5A–C; data Dryad doi:10.5061/dryad.k98sf7m8x). The prompt stated the biological
question, background subtraction, DAPI-defined nuclear regions and Costes thresholds, but did not name a plugin.
**a**, Pearson's R between Su(Hw) and γH2Av above the Costes thresholds, per image field, as published (values from the
authors' Figure 5 data sheet) and as measured by Agentic-J, which chose Fiji's Coloc 2 plugin. Bars, mean ± s.d.;
Δ, difference of means; P, exact two-sided permutation test. **b**, Agentic-J field values against published field values,
both sorted within each condition (the published sheet does not identify which value belongs to which image); dashed line,
identity. **c**, Increase in mean R under osmotic stress (NaCl − isotonic; 95% bootstrap CI, 10,000 resamples) for the
published analysis and for Agentic-J runs that differ in prompt, interaction mode and tool. Filled symbols, Pearson's R
above the Costes thresholds (the published metric); open symbols, Pearson's R over the whole nuclear region without
threshold, which runs given the original prompt (v1) reported because it did not specify the metric. Dotted line, published effect. All runs reproduce the increase;
the published effect size is matched when the published metric is used.
