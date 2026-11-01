#!/bin/sh
#![allow(unused)] /*
exec cargo +nightly -Zscript "$0" "$@"
*/

//! ```cargo
//! [package]
//! edition = "2024"
//! ```

use std::fs::File;
use std::io::{self, BufRead, Write};

fn main() -> io::Result<()> {
    println!("Extracting relativistic fat bands...");

    let proj_file = File::open("NiO.projwfc_bands.out").expect("Missing projwfc output!");
    let mut weights_ni = Vec::new();
    let mut weights_o = Vec::new();

    let (mut current_k, mut current_o_k) = (Vec::new(), Vec::new());
    let (mut cur_ni, mut cur_o) = (0.0, 0.0);
    let mut in_band = false;

    for line in io::BufReader::new(proj_file).lines().filter_map(Result::ok) {
        if line.starts_with(" k =") {
            if in_band {
                current_k.push(cur_ni);
                current_o_k.push(cur_o);
                cur_ni = 0.0;
                cur_o = 0.0;
            }
            if !current_k.is_empty() {
                weights_ni.push(current_k.clone());
                weights_o.push(current_o_k.clone());
                current_k.clear();
                current_o_k.clear();
            }
            in_band = false;
        } else if line.contains("==== e(") {
            if in_band {
                current_k.push(cur_ni);
                current_o_k.push(cur_o);
            }
            cur_ni = 0.0;
            cur_o = 0.0;
            in_band = true;
        } else if line.contains("* [#") {
            let segments: Vec<&str> = line.split("* [#").collect();
            for i in 1..segments.len() {
                let w: f64 = segments[i - 1]
                    .split_whitespace()
                    .last()
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0.0);
                let info = segments[i];
                if (info.contains("(Ni1)") || info.contains("(Ni2)")) && info.contains("(d_j") {
                    cur_ni += w;
                } else if info.contains("(O)") && info.contains("(p_j") {
                    cur_o += w;
                }
            }
        }
    }
    if in_band {
        current_k.push(cur_ni);
        current_o_k.push(cur_o);
        weights_ni.push(current_k);
        weights_o.push(current_o_k);
    }

    let bands_file = File::open("NiO_bands.dat.gnu").expect("Missing NiO_bands.dat.gnu!");
    let mut out = File::create("fatbands.dat")?;
    let (mut k_idx, mut band_idx) = (0, 0);

    for line in io::BufReader::new(bands_file)
        .lines()
        .filter_map(Result::ok)
    {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            if k_idx > 0 {
                band_idx += 1;
                k_idx = 0;
            }
            writeln!(out, "")?;
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            let (x, e) = (parts[0], parts[1]);
            let ni_w = weights_ni
                .get(k_idx)
                .and_then(|b| b.get(band_idx))
                .unwrap_or(&0.0);
            let o_w = weights_o
                .get(k_idx)
                .and_then(|b| b.get(band_idx))
                .unwrap_or(&0.0);
            writeln!(out, "{} {} {:.4} {:.4}", x, e, ni_w, o_w)?;
            k_idx += 1;
        }
    }
    Ok(())
}
