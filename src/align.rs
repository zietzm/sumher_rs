use std::collections::HashMap;

use crate::io::gwas::{GwasLine, GwasSumstats};
use crate::io::tagging::{TagInfo, TagLine};

pub fn check_aligned(gwas: &GwasSumstats, taginfo: &TagInfo) -> bool {
    gwas.predictors == taginfo.predictors
}

pub fn align_gwas_taginfo(gwas: &GwasSumstats, taginfo: &TagInfo) -> (GwasSumstats, TagInfo) {
    let gwas_predictor_to_idx: HashMap<String, usize> = gwas
        .results
        .iter()
        .enumerate()
        .map(|(i, x)| (x.predictor.clone(), i))
        .collect();

    let tag_predictor_to_idx: HashMap<String, usize> = taginfo
        .taglines
        .iter()
        .enumerate()
        .map(|(i, x)| (x.predictor.clone(), i))
        .collect();

    let common_predictors: Vec<String> = taginfo
        .predictors
        .clone()
        .into_iter()
        .filter(|x| gwas_predictor_to_idx.contains_key(x))
        .collect();

    let gwas_results: Vec<GwasLine> = common_predictors
        .iter()
        .map(|x| gwas_predictor_to_idx[x])
        .map(|x| gwas.results[x].clone())
        .collect();

    let taginfo_taglines: Vec<TagLine> = common_predictors
        .iter()
        .map(|x| tag_predictor_to_idx[x])
        .map(|x| taginfo.taglines[x].clone())
        .collect();

    (
        GwasSumstats::new(gwas.phenotype.clone(), gwas_results),
        TagInfo::new(
            taginfo_taglines,
            taginfo.annotation_names.clone(),
            taginfo.ssums.clone(),
        ),
    )
}

pub fn align_gwas_gwas(
    left: &GwasSumstats,
    right: &GwasSumstats,
    taginfo: &TagInfo,
) -> (GwasSumstats, GwasSumstats, TagInfo) {
    let left_predictor_to_idx: HashMap<String, usize> = left
        .results
        .iter()
        .enumerate()
        .map(|(i, x)| (x.predictor.clone(), i))
        .collect();

    let right_predictor_to_idx: HashMap<String, usize> = right
        .results
        .iter()
        .enumerate()
        .map(|(i, x)| (x.predictor.clone(), i))
        .collect();

    let tag_predictor_to_idx: HashMap<String, usize> = taginfo
        .taglines
        .iter()
        .enumerate()
        .map(|(i, x)| (x.predictor.clone(), i))
        .collect();

    let common_predictors: Vec<String> = left
        .predictors
        .clone()
        .into_iter()
        .filter(|x| right_predictor_to_idx.contains_key(x) && tag_predictor_to_idx.contains_key(x))
        .collect();

    let left_results: Vec<GwasLine> = common_predictors
        .iter()
        .map(|x| left_predictor_to_idx[x])
        .map(|x| left.results[x].clone())
        .collect();

    let right_results: Vec<GwasLine> = common_predictors
        .iter()
        .map(|x| right_predictor_to_idx[x])
        .map(|x| right.results[x].clone())
        .collect();

    let taginfo_taglines: Vec<TagLine> = common_predictors
        .iter()
        .map(|x| tag_predictor_to_idx[x])
        .map(|x| taginfo.taglines[x].clone())
        .collect();

    (
        GwasSumstats::new(left.phenotype.clone(), left_results),
        GwasSumstats::new(right.phenotype.clone(), right_results),
        TagInfo::new(
            taginfo_taglines,
            taginfo.annotation_names.clone(),
            taginfo.ssums.clone(),
        ),
    )
}
