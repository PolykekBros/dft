#!/bin/bash
#![allow(unused)] /*
exec cargo +nightly -Zscript "$0" "$@"
*/
use std::fs::File;
use std::io::{self, BufRead, Write};

// Clebsch-Gordan spin-up squared weights for l=2
// W_UP = (m_j + 2.5) / 5.0 for j=5/2 (m_j = -2.5, -1.5, -0.5, 0.5, 1.5, 2.5)
const W_UP_5_2: [f64; 6] = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0];

// W_UP = (2.5 - m_j) / 5.0 for j=3/2 (m_j = -1.5, -0.5, 0.5, 1.5)
const W_UP_3_2: [f64; 4] = [0.8, 0.6, 0.4, 0.2];

fn process_atom(atom_name: &str, file_j52: &str, file_j32: &str) -> io::Result<()> {
    let f52 = io::BufReader::new(File::open(file_j52)?);
    let f32 = io::BufReader::new(File::open(file_j32)?);

    let mut out = File::create(format!("{}_Sz_projected.dat", atom_name))?;
    writeln!(out, "# Energy(eV)  DOS_up  DOS_down  Net_Sz")?;

    // Read both j=5/2 and j=3/2 files line by line simultaneously
    for (line52, line32) in f52.lines().zip(f32.lines()) {
        let l52 = line52?;
        let l32 = line32?;
        // Skip header lines
        if l52.starts_with('#') || l32.starts_with('#') {
            continue;
        }

        let cols52: Vec<f64> = l52
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        let cols32: Vec<f64> = l32
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        // Col 0: Energy, Col 1: Total LDOS, Cols 2+: mj components
        if cols52.len() < 8 || cols32.len() < 6 {
            continue;
        }

        let energy = cols52[0];
        let mut dos_up = 0.0;
        let mut dos_down = 0.0;

        // Apply Clebsch-Gordan weights to j=5/2 mj components
        for i in 0..6 {
            let pdos = cols52[i + 2];
            dos_up += W_UP_5_2[i] * pdos;
            dos_down += (1.0 - W_UP_5_2[i]) * pdos;
        }

        // Apply Clebsch-Gordan weights to j=3/2 mj components
        for i in 0..4 {
            let pdos = cols32[i + 2];
            dos_up += W_UP_3_2[i] * pdos;
            dos_down += (1.0 - W_UP_3_2[i]) * pdos;
        }

        writeln!(
            out,
            "{:.4} {:.6} {:.6} {:.6}",
            energy,
            dos_up,
            dos_down,
            dos_up - dos_down
        )?;
    }
    Ok(())
}

fn main() {
    println!("Extracting exact Sz projections... don't mess up your filenames!");

    // Process Ni1
    if let Err(e) = process_atom(
        "Ni1",
        "NiO.pdos_atm#1(Ni1)_wfc#4(d_j2.5)",
        "NiO.pdos_atm#1(Ni1)_wfc#5(d_j1.5)",
    ) {
        eprintln!("Failed on Ni1: {}", e);
    }

    // Process Ni2
    if let Err(e) = process_atom(
        "Ni2",
        "NiO.pdos_atm#2(Ni2)_wfc#4(d_j2.5)",
        "NiO.pdos_atm#2(Ni2)_wfc#5(d_j1.5)",
    ) {
        eprintln!("Failed on Ni2: {}", e);
    }

    println!("Finished! Plot column 4 of the new .dat files to see the true magnetic splitting.");
}
