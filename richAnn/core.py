"""Core classes for richAnn package"""

import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)


class EnrichResult:
    """
    Class to store and manipulate enrichment analysis results
    
    Attributes:
    -----------
    result : pd.DataFrame
        Enrichment results table
    enrichment_type : str
        Type of enrichment analysis
    parameters : dict
        Analysis parameters
    """
    
    def __init__(self, 
                 result_df: pd.DataFrame, 
                 enrichment_type: str = "ORA",
                 parameters: dict = None):
        """Initialize EnrichResult object"""
        self.result = result_df.copy()
        self.enrichment_type = enrichment_type
        self.parameters = parameters or {}
        
        if self.result.index.name != 'Annot':
            self.result = self.result.reset_index(drop=True)
    
    def filter(self, 
               pvalue: float = None, 
               padj: float = None,
               min_count: int = 0,
               max_count: int = None,
               min_richfactor: float = 0.0,
               max_richfactor: float = None):
        """Filter enrichment results by various criteria"""
        df = self.result.copy()
        
        if pvalue is not None:
            df = df[df['Pvalue'] <= pvalue]
        if padj is not None:
            df = df[df['Padj'] <= padj]
        if min_count > 0:
            df = df[df['Count'] >= min_count]
        if max_count is not None:
            df = df[df['Count'] <= max_count]
        if min_richfactor > 0:
            df = df[df['RichFactor'] >= min_richfactor]
        if max_richfactor is not None:
            df = df[df['RichFactor'] <= max_richfactor]
        
        if len(df) == 0:
            logger.warning("No terms remaining after filtering")
        
        return EnrichResult(df, self.enrichment_type, self.parameters)
    
    def top(self, n: int = 10, orderby: str = 'Padj'):
        """Get top n enriched terms"""
        if orderby not in self.result.columns:
            raise ValueError(f"Column '{orderby}' not found in results")
        
        ascending = orderby in ['Pvalue', 'Padj']
        df = self.result.sort_values(orderby, ascending=ascending)
        
        return EnrichResult(df.head(n), self.enrichment_type, self.parameters)
    
    def detail(self):
        """Extract detailed gene-term relationships"""
        records = []
        for _, row in self.result.iterrows():
            genes = str(row['GeneID']).replace(',', ';').split(';')
            for gene in genes:
                gene = gene.strip()
                if gene:
                    records.append({
                        'Gene': gene,
                        'Annot': row['Annot'],
                        'Term': row['Term'],
                        'Pvalue': row['Pvalue'],
                        'Padj': row['Padj'],
                        'RichFactor': row['RichFactor'],
                        'Count': row['Count']
                    })
        
        return pd.DataFrame(records)
    
    def summary(self):
        """Print comprehensive summary statistics"""
        df = self.result

        print(f"\n{'='*70}")
        print(f"  richAnn Enrichment Result ({self.enrichment_type})")
        print(f"{'='*70}")
        print(f"\nDataset Overview:")
        print(f"  Total terms: {len(df)}")
        print(f"  Significant (Padj < 0.05): {len(df[df['Padj'] < 0.05])}")
        print(f"  Significant (Padj < 0.01): {len(df[df['Padj'] < 0.01])}")

        # GSEA-specific statistics
        if self.enrichment_type == "GSEA":
            if 'ES' in df.columns:
                print(f"\nEnrichment Score Statistics:")
                print(f"  ES range: {df['ES'].min():.3f} to {df['ES'].max():.3f}")
                print(f"  Mean ES: {df['ES'].mean():.3f}")

            if 'NES' in df.columns:
                print(f"\nNormalized Enrichment Score Statistics:")
                print(f"  NES range: {df['NES'].min():.3f} to {df['NES'].max():.3f}")
                print(f"  Mean NES: {df['NES'].mean():.3f}")
                print(f"  Highly enriched (|NES| > 2): {len(df[df['NES'].abs() > 2])}")

            print(f"\nGene Set Size Statistics:")
            print(f"  Range: {df['Count'].min()} - {df['Count'].max()}")
            print(f"  Mean: {df['Count'].mean():.1f}")
            print(f"  Median: {df['Count'].median():.1f}")

            # Top terms by NES
            print(f"\nTop 5 Gene Sets by |NES|:")
            df_abs_nes = df.copy()
            df_abs_nes['abs_NES'] = df_abs_nes['NES'].abs()
            top5 = df_abs_nes.nlargest(5, 'abs_NES')[['Term', 'NES', 'Pvalue', 'Padj']]
            print(top5.to_string(index=False))

            print(f"\nTop 5 Gene Sets by Significance:")
            top5_sig = df.nsmallest(5, 'Padj')[['Term', 'NES', 'Pvalue', 'Padj']]
            print(top5_sig.to_string(index=False))

        # ORA/GO/KEGG-specific statistics
        else:
            if 'RichFactor' in df.columns:
                print(f"  Highly enriched (RichFactor > 2): {len(df[df['RichFactor'] > 2])}")

                print(f"\nEnrichment Statistics:")
                print(f"  RichFactor range: {df['RichFactor'].min():.3f} - {df['RichFactor'].max():.3f}")
                print(f"  Mean RichFactor: {df['RichFactor'].mean():.3f}")
                print(f"  Median RichFactor: {df['RichFactor'].median():.3f}")

            print(f"\nGene Count Statistics:")
            print(f"  Range: {df['Count'].min()} - {df['Count'].max()}")
            print(f"  Mean: {df['Count'].mean():.1f}")
            print(f"  Median: {df['Count'].median():.1f}")

            if 'OddsRatio' in df.columns:
                print(f"\nOdds Ratio Statistics:")
                print(f"  Range: {df['OddsRatio'].min():.3f} - {df['OddsRatio'].max():.3f}")
                print(f"  Mean: {df['OddsRatio'].mean():.3f}")

            if 'RichFactor' in df.columns:
                print(f"\nTop 5 Terms by RichFactor:")
                top5 = df.nlargest(5, 'RichFactor')[['Term', 'RichFactor', 'Count', 'Padj']]
                print(top5.to_string(index=False))

                print(f"\nTop 5 Terms by Significance:")
                top5_sig = df.nsmallest(5, 'Padj')[['Term', 'RichFactor', 'Count', 'Padj']]
                print(top5_sig.to_string(index=False))

        # Parameters (common for all types)
        if self.parameters:
            print(f"\nAnalysis Parameters:")
            for key, value in self.parameters.items():
                print(f"  {key}: {value}")

        print(f"{'='*70}\n")
    
    def to_csv(self, filename: str, **kwargs):
        """Export results to CSV file"""
        self.result.to_csv(filename, index=False, **kwargs)
        logger.info(f"Results exported to {filename}")
    
    def to_excel(self, filename: str, sheet_name: str = 'Enrichment', **kwargs):
        """Export results to Excel file"""
        self.result.to_excel(filename, sheet_name=sheet_name, index=False, **kwargs)
        logger.info(f"Results exported to {filename}")
    
    def __len__(self):
        return len(self.result)
    
    def __repr__(self):
        n_terms = len(self.result)
        n_sig = len(self.result[self.result['Padj'] < 0.05])
        return (f"EnrichResult ({self.enrichment_type})\n"
                f"{n_terms} terms, {n_sig} significant (Padj < 0.05)\n"
                f"Use .summary() for detailed statistics")
    
    def __getitem__(self, key):
        """Allow indexing like a DataFrame"""
        return self.result[key]
