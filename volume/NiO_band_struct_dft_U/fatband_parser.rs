#!/bin/sh
#![allow(unused)] /*
exec cargo +nightly -Zscript "$0" "$@"
*/

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, Write};

fn main() -> io::Result<()> {
    println!("Extracting relativistic fat bands... using bulletproof parsing!");
    
    let proj_file = match File::open("NiO.projwfc_bands.out") {
        Ok(f) => f,
        Err(_) => {
            eprintln!("Error: NiO.projwfc_bands.out not found!");
            return Ok(());
        }
    };
    
    let mut ni_states = HashSet::new();
    let mut o_states = HashSet::new();
    
    let mut weights_ni = Vec::new();
    let mut weights_o = Vec::new();
    
    let (mut current_k, mut current_o_k) = (Vec::new(), Vec::new());
    let (mut cur_ni, mut cur_o) = (0.0, 0.0);
    let mut in_band = false;
    
    for line_result in io::BufReader::new(proj_file).lines() {
        let line = line_result?;
        
        // 1. Bulletproof state mapping
        if line.contains("state #") && line.contains("atom") {
            let parts: Vec<&str> = line.split("state #").collect();
            if parts.len() > 1 {
                let state_part = parts[1].split(':').next().unwrap_or("").trim();
                if let Ok(state_num) = state_part.parse::<usize>() {
                    // Match any Ni d-orbital and O p-orbital regardless of specific naming
                    if line.contains("Ni") && (line.contains("l=2") || line.contains("d")) {
                        ni_states.insert(state_num);
                    } else if line.contains("O") && (line.contains("l=1") || line.contains("p")) {
                        o_states.insert(state_num);
                    }
                }
            }
            continue;
        }

        // 2. Process k-points and energy blocks
        if line.starts_with(" k =") {
            if in_band {
                current_k.push(cur_ni);
                current_o_k.push(cur_o);
                weights_ni.push(current_k.clone());
                weights_o.push(current_o_k.clone());
                current_k.clear(); 
                current_o_k.clear();
                cur_ni = 0.0; 
                cur_o = 0.0;
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
        } 
        // 3. Bulletproof weight extraction based on bracket delimiters
        else if line.contains('[') && (line.contains("*") || line.contains("+")) {
            let segments: Vec<&str> = line.split('[').collect();
            for i in 1..segments.len() {
                let inside = segments[i].split(']').next().unwrap_or("");
                let state_str = inside.replace("#", "").trim().to_string();
                
                if let Ok(state_num) = state_str.parse::<usize>() {
                    let w_str = segments[i-1].replace("*", " ").replace("+", " ").replace("=", " ").replace("|psi|^2", " ");
                    let w: f64 = w_str.split_whitespace().last().unwrap_or("0").parse().unwrap_or(0.0);
                    
                    if ni_states.contains(&state_num) { cur_ni += w; }
                    if o_states.contains(&state_num) { cur_o += w; }
                }
            }
        }
    }
    
    // Push the very last band block
    if in_band {
        current_k.push(cur_ni); 
        current_o_k.push(cur_o);
        weights_ni.push(current_k); 
        weights_o.push(current_o_k);
    }
    
    // Diagnostics to ensure data was actually captured
    println!("Diagnostic: Found {} Ni states and {} O states.", ni_states.len(), o_states.len());
    println!("Diagnostic: Parsed {} k-points.", weights_ni.len());
    if let Some(first) = weights_ni.first() {
        println!("Diagnostic: First k-point has {} bands.", first.len());
    }
    
    let bands_file = match File::open("NiO_bands.dat.gnu") {
        Ok(f) => f,
        Err(_) => {
            eprintln!("Error: NiO_bands.dat.gnu not found!");
            return Ok(());
        }
    };
    
    let mut out = File::create("fatbands.dat")?;
    let (mut k_idx, mut band_idx) = (0, 0);
    let mut points_written = 0;
    
    // 4. Align with Gnuplot structure
    for line in io::BufReader::new(bands_file).lines().filter_map(Result::ok) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            if k_idx > 0 { band_idx += 1; k_idx = 0; }
            writeln!(out, "")?;
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            let (x, e) = (parts[0], parts[1]);
            let ni_w = weights_ni.get(k_idx).and_then(|b| b.get(band_idx)).unwrap_or(&0.0);
            let o_w = weights_o.get(k_idx).and_then(|b| b.get(band_idx)).unwrap_or(&0.0);
            writeln!(out, "{} {} {:.4} {:.4}", x, e, ni_w, o_w)?;
            
            if *ni_w > 0.01 || *o_w > 0.01 { points_written += 1; }
            k_idx += 1;
        }
    }
    
    println!("Diagnostic: Wrote fatbands.dat with {} non-zero orbital points.", points_written);
    Ok(())
}
