"""
Tests for enrichment analysis functions
"""

import pytest
import pandas as pd
import numpy as np
from richAnn import richGO, richKEGG, richGSEA
from richAnn.core import EnrichResult


class TestRichGO:
    """Tests for GO enrichment analysis"""

    def test_basic_enrichment(self, sample_genes, sample_go_data):
        """Test basic GO enrichment"""
        result = richGO(sample_genes, sample_go_data, ontology="BP", pvalue=1.0)

        assert isinstance(result, EnrichResult)
        assert len(result) > 0
        assert 'Annot' in result.result.columns
        assert 'Term' in result.result.columns
        assert 'Pvalue' in result.result.columns
        assert 'Padj' in result.result.columns

    def test_ontology_filter(self, sample_genes, sample_go_data):
        """Test ontology filtering"""
        result = richGO(sample_genes, sample_go_data, ontology="BP", pvalue=1.0)

        # All results should be BP
        if len(result) > 0:
            assert all(result.result['Annot'].str.startswith('GO:'))

    def test_invalid_ontology(self, sample_genes, sample_go_data):
        """Test invalid ontology raises error"""
        with pytest.raises(ValueError, match="Ontology must be"):
            richGO(sample_genes, sample_go_data, ontology="INVALID")

    def test_empty_genes(self, sample_go_data):
        """Test empty gene list raises error"""
        with pytest.raises(ValueError, match="No valid genes"):
            richGO([], sample_go_data, ontology="BP")

    def test_invalid_pvalue(self, sample_genes, sample_go_data):
        """Test invalid p-value raises error"""
        with pytest.raises(ValueError, match="pvalue must be between"):
            richGO(sample_genes, sample_go_data, ontology="BP", pvalue=1.5)

    def test_case_insensitive(self, sample_go_data):
        """Test case-insensitive gene matching"""
        genes_lower = ["brca1", "brca2", "tp53"]
        result = richGO(genes_lower, sample_go_data, ontology="BP",
                       pvalue=1.0, case_sensitive=False)

        assert isinstance(result, EnrichResult)

    def test_case_sensitive(self, sample_go_data):
        """Test case-sensitive gene matching"""
        genes_lower = ["brca1", "brca2", "tp53"]

        # Should raise error if no matches with case sensitivity
        with pytest.raises(ValueError, match="No query genes found"):
            richGO(genes_lower, sample_go_data, ontology="BP",
                  pvalue=1.0, case_sensitive=True)


class TestRichKEGG:
    """Tests for KEGG pathway enrichment"""

    def test_basic_enrichment(self, sample_genes, sample_kegg_data):
        """Test basic KEGG enrichment"""
        result = richKEGG(sample_genes, sample_kegg_data, pvalue=1.0)

        assert isinstance(result, EnrichResult)
        assert 'Pathway' in result.result['Annot'].iloc[0] or 'hsa' in result.result['Annot'].iloc[0]

    def test_empty_genes(self, sample_kegg_data):
        """Test empty gene list raises error"""
        with pytest.raises(ValueError, match="No valid genes"):
            richKEGG([], sample_kegg_data)


class TestRichGSEA:
    """Tests for GSEA"""

    def test_basic_gsea(self, sample_gene_scores, sample_geneset_db):
        """Test basic GSEA"""
        result = richGSEA(sample_gene_scores, sample_geneset_db,
                         nperm=100, min_size=2, max_size=100)

        assert isinstance(result, EnrichResult)
        assert len(result) > 0
        assert 'ES' in result.result.columns
        assert 'NES' in result.result.columns
        assert 'LeadingEdge' in result.result.columns

    def test_empty_scores(self, sample_geneset_db):
        """Test empty gene scores raises error"""
        with pytest.raises(ValueError, match="cannot be empty"):
            richGSEA({}, sample_geneset_db)

    def test_invalid_nperm(self, sample_gene_scores, sample_geneset_db):
        """Test invalid nperm raises error"""
        with pytest.raises(ValueError, match="must be positive"):
            richGSEA(sample_gene_scores, sample_geneset_db, nperm=0)

    def test_weighted_vs_unweighted(self, sample_gene_scores, sample_geneset_db):
        """Test weighted vs unweighted GSEA"""
        result_weighted = richGSEA(sample_gene_scores, sample_geneset_db,
                                   nperm=50, weighted_score=True, random_state=42)
        result_unweighted = richGSEA(sample_gene_scores, sample_geneset_db,
                                     nperm=50, weighted_score=False, random_state=42)

        assert isinstance(result_weighted, EnrichResult)
        assert isinstance(result_unweighted, EnrichResult)
        # Scores should differ between weighted and unweighted
        if len(result_weighted) > 0 and len(result_unweighted) > 0:
            assert not np.allclose(result_weighted.result['ES'].values,
                                  result_unweighted.result['ES'].values)


class TestEnrichResultMethods:
    """Tests for EnrichResult methods"""

    def test_filter_by_pvalue(self, sample_enrich_result):
        """Test filtering by p-value"""
        filtered = sample_enrich_result.filter(pvalue=0.05)

        assert isinstance(filtered, EnrichResult)
        if len(filtered) > 0:
            assert all(filtered.result['Pvalue'] <= 0.05)

    def test_filter_by_count(self, sample_enrich_result):
        """Test filtering by gene count"""
        filtered = sample_enrich_result.filter(min_count=3)

        assert isinstance(filtered, EnrichResult)
        if len(filtered) > 0:
            assert all(filtered.result['Count'] >= 3)

    def test_top_n(self, sample_enrich_result):
        """Test getting top N results"""
        top5 = sample_enrich_result.top(n=5)

        assert isinstance(top5, EnrichResult)
        assert len(top5) <= 5

    def test_detail(self, sample_enrich_result):
        """Test detail method"""
        detail_df = sample_enrich_result.detail()

        assert isinstance(detail_df, pd.DataFrame)
        assert 'Gene' in detail_df.columns
        assert 'Annot' in detail_df.columns
        assert 'Term' in detail_df.columns

    def test_summary(self, sample_enrich_result, capsys):
        """Test summary method prints output"""
        sample_enrich_result.summary()

        captured = capsys.readouterr()
        assert 'richAnn Enrichment Result' in captured.out
        assert 'Total terms:' in captured.out
