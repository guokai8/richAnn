"""
Tests for data loading utilities
"""

import pytest
import pandas as pd
import tempfile
import os
from pathlib import Path
from richAnn.data import (
    load_gmt, create_custom_annotation, validate_annotation_format
)


class TestLoadGMT:
    """Tests for GMT file loading"""

    def test_load_valid_gmt(self):
        """Test loading a valid GMT file"""
        # Create a temporary GMT file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gmt', delete=False) as f:
            f.write("GENESET1\tDescription1\tGENE1\tGENE2\tGENE3\n")
            f.write("GENESET2\tDescription2\tGENE4\tGENE5\n")
            f.write("GENESET3\tDescription3\tGENE1\tGENE6\tGENE7\tGENE8\n")
            temp_path = f.name

        try:
            df = load_gmt(temp_path)

            assert isinstance(df, pd.DataFrame)
            assert 'GeneSet' in df.columns
            assert 'GeneSetName' in df.columns
            assert 'GeneID' in df.columns
            assert df['GeneSet'].nunique() == 3
            assert 'GENE1' in df['GeneID'].values
            assert 'GENESET1' in df['GeneSet'].values

        finally:
            os.unlink(temp_path)

    def test_load_empty_gmt(self):
        """Test loading an empty GMT file raises error"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gmt', delete=False) as f:
            temp_path = f.name

        try:
            with pytest.raises(ValueError, match="No valid gene sets found"):
                load_gmt(temp_path)
        finally:
            os.unlink(temp_path)

    def test_load_nonexistent_file(self):
        """Test loading non-existent file raises error"""
        with pytest.raises(FileNotFoundError):
            load_gmt("nonexistent_file.gmt")

    def test_custom_column_names(self):
        """Test using custom column names"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gmt', delete=False) as f:
            f.write("SET1\tDesc1\tG1\tG2\n")
            temp_path = f.name

        try:
            df = load_gmt(temp_path, annot_col="Pathway", name_col="Name", gene_col="Gene")

            assert 'Pathway' in df.columns
            assert 'Name' in df.columns
            assert 'Gene' in df.columns

        finally:
            os.unlink(temp_path)


class TestCreateCustomAnnotation:
    """Tests for creating custom annotations"""

    def test_create_basic_annotation(self):
        """Test creating annotation from dictionary"""
        gene_sets = {
            "SET1": ["GENE1", "GENE2", "GENE3"],
            "SET2": ["GENE4", "GENE5"],
            "SET3": ["GENE1", "GENE6", "GENE7"]
        }

        df = create_custom_annotation(gene_sets)

        assert isinstance(df, pd.DataFrame)
        assert 'Annot' in df.columns
        assert 'Term' in df.columns
        assert 'GeneID' in df.columns
        assert df['Annot'].nunique() == 3
        assert 'GENE1' in df['GeneID'].values

    def test_empty_gene_sets(self):
        """Test empty gene sets dictionary raises error"""
        with pytest.raises(ValueError, match="cannot be empty"):
            create_custom_annotation({})

    def test_skip_empty_sets(self):
        """Test that empty gene sets are skipped"""
        gene_sets = {
            "SET1": ["GENE1", "GENE2"],
            "SET2": [],  # Empty
            "SET3": ["GENE3"]
        }

        df = create_custom_annotation(gene_sets)

        # SET2 should be skipped
        assert df['Annot'].nunique() == 2
        assert 'SET2' not in df['Annot'].values

    def test_custom_column_names(self):
        """Test using custom column names"""
        gene_sets = {"SET1": ["G1", "G2"]}

        df = create_custom_annotation(
            gene_sets,
            annot_col="Pathway",
            name_col="PathName",
            gene_col="Gene"
        )

        assert 'Pathway' in df.columns
        assert 'PathName' in df.columns
        assert 'Gene' in df.columns


class TestValidateAnnotation:
    """Tests for annotation validation"""

    def test_validate_valid_annotation(self):
        """Test validation of valid annotation"""
        df = pd.DataFrame({
            'GeneID': ['G1', 'G2', 'G3'],
            'Annot': ['A1', 'A1', 'A2'],
            'Term': ['Term1', 'Term1', 'Term2']
        })

        result = validate_annotation_format(df, ['GeneID', 'Annot', 'Term'])

        assert result['valid'] == True
        assert result['n_rows'] == 3
        assert result['n_annotations'] == 2
        assert result['n_genes'] == 3

    def test_validate_missing_columns(self):
        """Test validation fails with missing columns"""
        df = pd.DataFrame({
            'GeneID': ['G1', 'G2'],
            'Annot': ['A1', 'A1']
        })

        result = validate_annotation_format(df, ['GeneID', 'Annot', 'Term'])

        assert result['valid'] == False
        assert len(result['errors']) > 0
        assert 'Missing required columns' in result['errors'][0]

    def test_validate_empty_dataframe(self):
        """Test validation fails with empty DataFrame"""
        df = pd.DataFrame(columns=['GeneID', 'Annot', 'Term'])

        result = validate_annotation_format(df, ['GeneID', 'Annot', 'Term'])

        assert result['valid'] == False
        assert 'empty' in result['errors'][0].lower()

    def test_validate_with_nulls(self):
        """Test validation warns about null values"""
        df = pd.DataFrame({
            'GeneID': ['G1', None, 'G3'],
            'Annot': ['A1', 'A1', 'A2'],
            'Term': ['Term1', 'Term1', 'Term2']
        })

        result = validate_annotation_format(df, ['GeneID', 'Annot', 'Term'])

        assert len(result['warnings']) > 0
        assert any('null' in w.lower() for w in result['warnings'])

    def test_validate_with_duplicates(self):
        """Test validation warns about duplicates"""
        df = pd.DataFrame({
            'GeneID': ['G1', 'G1', 'G2'],
            'Annot': ['A1', 'A1', 'A1'],
            'Term': ['Term1', 'Term1', 'Term1']
        })

        result = validate_annotation_format(df, ['GeneID', 'Annot', 'Term'])

        assert len(result['warnings']) > 0
        assert any('duplicate' in w.lower() for w in result['warnings'])
