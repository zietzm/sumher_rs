use anyhow::{Context, Result};
use rusqlite::{Connection, Transaction};
use serde_rusqlite::to_params_named;
use std::collections::HashSet;
use std::path::Path;

use crate::hsq::{HsqResult, RgResult};

pub struct DbConnection {
    pub conn: Connection,
}

pub struct Progress {
    pub finished_h2: HashSet<String>,
    pub finished_rg: HashSet<(String, String)>,
}

impl DbConnection {
    pub fn new<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let preexisting = path.as_ref().is_file();

        let conn = Connection::open(path)?;

        if !preexisting {
            conn.execute(
                "CREATE TABLE IF NOT EXISTS h2 (
                    phenotype TEXT NOT NULL,
                    component TEXT NOT NULL,
                    h2 FLOAT NOT NULL,
                    h2_se FLOAT NOT NULL
                )",
                (),
            )?;

            conn.execute(
                "CREATE TABLE IF NOT EXISTS rg (
                    phenotype1 TEXT NOT NULL,
                    phenotype2 TEXT NOT NULL,
                    component TEXT NOT NULL,
                    rg FLOAT NOT NULL,
                    rg_se FLOAT NOT NULL
                )",
                (),
            )?;
        }

        Ok(Self { conn })
    }

    pub fn get_progress(&self) -> Result<Progress> {
        let finished_h2 = self
            .conn
            .prepare("SELECT phenotype FROM h2 WHERE component = 'Her_All'")?
            .query_map([], |row| row.get(0))?
            .map(|x| x.unwrap())
            .collect::<HashSet<_>>();

        let finished_rg = self
            .conn
            .prepare("SELECT phenotype1, phenotype2 FROM rg WHERE component = 'Cor_All'")?
            .query_map([], |row| Ok((row.get(0), row.get(1))))?
            .map(|x| x.unwrap())
            .map(|(x, y)| (x.unwrap(), y.unwrap()))
            .collect::<HashSet<_>>();

        Ok(Progress {
            finished_h2,
            finished_rg,
        })
    }

    pub fn write_h2(&mut self, rows: &[HsqResult]) -> Result<()> {
        let tx = self.conn.transaction()?;
        insert_hsq_data(&tx, rows)?;
        tx.commit()
            .context(format!("Failed to commit h2 data {:?}", rows[0]))?;

        Ok(())
    }

    pub fn write_rg(&mut self, rows: &[RgResult]) -> Result<()> {
        let tx = self.conn.transaction()?;
        insert_rg_data(&tx, rows)?;
        tx.commit()
            .context(format!("Failed to commit h2 data {:?}", rows[0]))?;

        Ok(())
    }
}

fn insert_hsq_data(tx: &Transaction, rows: &[HsqResult]) -> Result<()> {
    let mut stmt = tx.prepare_cached(
        "INSERT INTO h2 (phenotype, component, h2, h2_se) VALUES 
        (:phenotype, :component, :estimate, :se)",
    )?;

    for row in rows.iter() {
        stmt.execute(to_params_named(row)?.to_slice().as_slice())
            .context(format!("Failed to insert h2 data for {:?}", row))?;
    }

    Ok(())
}

fn insert_rg_data(tx: &Transaction, rows: &[RgResult]) -> Result<()> {
    let mut stmt = tx.prepare_cached(
        "INSERT INTO rg (phenotype1, phenotype2, component, rg, rg_se) VALUES 
        (:phenotype1, :phenotype2, :component, :estimate, :se)",
    )?;

    for row in rows.iter() {
        stmt.execute(to_params_named(row)?.to_slice().as_slice())
            .context(format!("Failed to insert rg data for {:?}", row))?;
    }

    Ok(())
}
