"""
Tests for utility functions
"""

import pytest
import pandas as pd
import numpy as np
from richAnn.utils import (
    _normalize_genes, _validate_annotation, _calculate_effect_sizes,
    _validate_pvalue, _validate_size_params, _validate_positive_int,
    compareResult
)
from richAnn.core import EnrichResult


class TestNormalizeGenes:
    """Tests for gene normalization"""

    def test_normalize_list(self):
        """Test normalizing a list of genes"""
        genes = ["BRCA1", "BRCA2", "TP53"]
        result = _normalize_genes(genes, case_sensitive=False)

        assert isinstance(result, set)
        assert "BRCA1" in result
        assert len(result) == 3

    def test_normalize_case_insensitive(self):
        """Test case-insensitive normalization"""
        genes = ["brca1", "BRCA2", "TpS3"]
        result = _normalize_genes(genes, case_sensitive=False)

        assert "BRCA1" in result
        assert "BRCA2" in result
        assert "TPS3" in result

    def test_normalize_case_sensitive(self):
        """Test case-sensitive normalization"""
        genes = ["brca1", "BRCA2"]
        result = _normalize_genes(genes, case_sensitive=True)

        assert "brca1" in result
        assert "BRCA2" in result
        assert "BRCA1" not in result

    def test_normalize_set(self):
        """Test normalizing a set"""
        genes = {"BRCA1", "BRCA2"}
        result = _normalize_genes(genes, case_sensitive=True)

        assert isinstance(result, set)
        assert len(result) == 2

    def test_normalize_array(self):
        """Test normalizing numpy array"""
        genes = np.array(["BRCA1", "BRCA2", "TP53"])
        result = _normalize_genes(genes, case_sensitive=True)

        assert isinstance(result, set)
        assert len(result) == 3

    def test_normalize_strips_whitespace(self):
        """Test that whitespace is stripped"""
        genes = [" BRCA1 ", "BRCA2\t", "\nTP53"]
        result = _normalize_genes(genes, case_sensitive=True)

        assert "BRCA1" in result
        assert " BRCA1 " not in result


class TestValidateAnnotation:
    """Tests for annotation validation"""

    def test_validate_valid_annotation(self):
        """Test validation passes for valid annotation"""
        df = pd.DataFrame({
            'GeneID': ['G1', 'G2'],
            'Annot': ['A1', 'A1'],
            'Term': ['T1', 'T1']
        })

        # Should not raise any exception
        _validate_annotation(df, ['GeneID', 'Annot', 'Term'])

    def test_validate_missing_columns(self):
        """Test validation fails with missing columns"""
        df = pd.DataFrame({
            'GeneID': ['G1', 'G2'],
            'Annot': ['A1', 'A1']
        })

        with pytest.raises(ValueError, match="missing required columns"):
            _validate_annotation(df, ['GeneID', 'Annot', 'Term'])

    def test_validate_not_dataframe(self):
        """Test validation fails if not DataFrame"""
        with pytest.raises(TypeError, match="must be a pandas DataFrame"):
            _validate_annotation([], ['GeneID'])

    def test_validate_empty_dataframe(self):
        """Test validation fails for empty DataFrame"""
        df = pd.DataFrame(columns=['GeneID', 'Annot'])

        with pytest.raises(ValueError, match="empty"):
            _validate_annotation(df, ['GeneID', 'Annot'])


class TestCalculateEffectSizes:
    """Tests for effect size calculations"""

    def test_basic_calculation(self):
        """Test basic effect size calculation"""
        result = _calculate_effect_sizes(k=5, n=10, M=20, N=100)

        assert 'RichFactor' in result
        assert 'FoldEnrich' in result
        assert 'OddsRatio' in result
        assert 'LogOddsRatio' in result

        # RichFactor = (k/n) / (M/N) = (5/10) / (20/100) = 0.5 / 0.2 = 2.5
        assert result['RichFactor'] == pytest.approx(2.5, rel=0.01)

    def test_no_division_by_zero(self):
        """Test that division by zero is handled"""
        result = _calculate_effect_sizes(k=0, n=10, M=20, N=100)

        assert result['RichFactor'] >= 0
        assert not np.isnan(result['RichFactor'])


class TestValidationFunctions:
    """Tests for validation helper functions"""

    def test_validate_pvalue_valid(self):
        """Test p-value validation with valid values"""
        _validate_pvalue(0.05)
        _validate_pvalue(0.001)
        _validate_pvalue(1.0)

    def test_validate_pvalue_invalid(self):
        """Test p-value validation with invalid values"""
        with pytest.raises(ValueError):
            _validate_pvalue(1.5)

        with pytest.raises(ValueError):
            _validate_pvalue(0.0)

        with pytest.raises(ValueError):
            _validate_pvalue(-0.1)

    def test_validate_size_params_valid(self):
        """Test size parameter validation with valid values"""
        _validate_size_params(5, 100)
        _validate_size_params(0, 500)

    def test_validate_size_params_invalid(self):
        """Test size parameter validation with invalid values"""
        with pytest.raises(ValueError):
            _validate_size_params(100, 50)  # max < min

        with pytest.raises(ValueError):
            _validate_size_params(-5, 100)  # negative min

    def test_validate_positive_int_valid(self):
        """Test positive integer validation with valid values"""
        _validate_positive_int(1, "test")
        _validate_positive_int(1000, "test")

    def test_validate_positive_int_invalid(self):
        """Test positive integer validation with invalid values"""
        with pytest.raises(ValueError):
            _validate_positive_int(0, "test")

        with pytest.raises(ValueError):
            _validate_positive_int(-5, "test")

        with pytest.raises(TypeError):
            _validate_positive_int(1.5, "test")


class TestCompareResult:
    """Tests for compareResult function"""

    def test_compare_basic(self):
        """Test basic result comparison"""
        result1 = EnrichResult(
            pd.DataFrame({
                'Annot': ['GO:001', 'GO:002'],
                'Term': ['Term1', 'Term2'],
                'Pvalue': [0.01, 0.02],
                'Padj': [0.05, 0.06],
                'Count': [5, 3],
                'RichFactor': [2.0, 1.5]
            }),
            enrichment_type="GO"
        )

        result2 = EnrichResult(
            pd.DataFrame({
                'Annot': ['GO:001', 'GO:003'],
                'Term': ['Term1', 'Term3'],
                'Pvalue': [0.005, 0.03],
                'Padj': [0.03, 0.07],
                'Count': [6, 4],
                'RichFactor': [2.5, 1.8]
            }),
            enrichment_type="GO"
        )

        result_dict = {'Sample1': result1, 'Sample2': result2}
        combined = compareResult(result_dict)

        assert isinstance(combined, pd.DataFrame)
        assert 'Sample' in combined.columns
        assert len(combined) == 4  # 2 + 2 results
        assert 'Sample1' in combined['Sample'].values
        assert 'Sample2' in combined['Sample'].values

    def test_compare_empty_dict(self):
        """Test comparing empty dictionary raises error"""
        with pytest.raises(ValueError, match="cannot be empty"):
            compareResult({})

    def test_compare_invalid_type(self):
        """Test comparing with non-EnrichResult raises error"""
        with pytest.raises(TypeError):
            compareResult({'Sample1': pd.DataFrame()})
