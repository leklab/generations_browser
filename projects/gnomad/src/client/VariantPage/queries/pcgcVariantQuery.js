export default `
query PcgcVariant($variantId: String!) {
  variant(variantId: $variantId) {
    alt
    chrom
    pos
    ref
    variantId
    xpos
    ... on VariantDetails {
      exome {
        ac
        an
        ac_hom
        populations {
          id
          ac
          an
          ac_hom
        }
      }
      sortedTranscriptConsequences {
        canonical
        gene_id
        gene_symbol
        hgvs
        hgvsc
        hgvsp
        lof
        lof_flags
        lof_filter
        lof_info
        major_consequence
        polyphen_prediction
        sift_prediction
        transcript_id
      }
    }
  }
}
`

/*
export default `
query PcgcVariant($variantId: String!, $datasetId: DatasetsSupportingFetchVariantDetails!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    alt
    chrom
    pos
    ref
    variantId
    xpos
    ... on GnomadVariantDetails {
      exome {
        ac
        an
        ac_hom
        populations {
          id
          ac
          an
        }

      }
      sortedTranscriptConsequences {
        canonical
        gene_id
        gene_symbol
        hgvs
        hgvsc
        hgvsp
        lof
        lof_flags
        lof_filter
        lof_info
        major_consequence
        polyphen_prediction
        sift_prediction
        transcript_id
      }
    }
  }
}
`
*/