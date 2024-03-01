// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Process the SUBSET_JSON option.

use enclone_core::defs::{EncloneControl, ExactClonotype};
use io_utils::{open_for_write_new, open_maybe_compressed, path_exists};
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;
use vector_utils::{bin_member, unique_sort};

#[derive(Serialize, Deserialize)]
struct AnnotationWithDataset<'a> {
    dataset: Option<&'a str>,
    #[serde(flatten)]
    data: ContigAnnotation,
}

pub fn subset_json(
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    exacts: &Vec<Vec<usize>>,
    ann: &str,
) -> Result<(), String> {
    if ctl.gen_opt.subset_json.is_empty() {
        return Ok(());
    }

    let mut barcode_li = Vec::<(&str, usize)>::new();
    for l in 0..exacts.len() {
        for u in 0..exacts[l].len() {
            let ex = &exact_clonotypes[exacts[l][u]];
            for j in 0..ex.clones.len() {
                barcode_li.push((
                    ex.clones[j][0].barcode.as_str(),
                    ex.clones[j][0].dataset_index,
                ));
            }
        }
    }
    unique_sort(&mut barcode_li);

    let annotations: Vec<_> =
        std::iter::zip(&ctl.origin_info.dataset_path, &ctl.origin_info.dataset_id)
            .enumerate()
            .flat_map(|(li, (ds_path, ds_id))| {
                let mut json_path = format!("{}/{}", ds_path, ann);
                if !path_exists(&json_path) {
                    json_path = format!("{}.lz4", json_path);
                }
                let mut contents = String::new();
                open_maybe_compressed(&json_path)
                    .read_to_string(&mut contents)
                    .unwrap();

                serde_json::Deserializer::from_str(&contents)
                    .into_iter::<AnnotationWithDataset>()
                    .map(Result::unwrap)
                    .filter(|ann| bin_member(&barcode_li, &(&ann.data.barcode, li)))
                    .map(|ann| AnnotationWithDataset {
                        dataset: Some(ds_id),
                        data: ann.data,
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

    serde_json::to_writer_pretty(open_for_write_new![&ctl.gen_opt.subset_json], &annotations)
        .map_err(|e| e.to_string())
}
