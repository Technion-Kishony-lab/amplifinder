"""Step 14: Export results to CSV."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.data_types import RecordTypedDF, AnalyzedTnJc2
from amplifinder.steps.base import Step
from amplifinder.logger import info


class ExportStep(Step[None]):
    """Export analyzed candidates to CSV files.
    
    Creates two CSV files:
    1. ISJC2.csv - All analyzed candidates
    2. candidate_amplifications.csv - Filtered candidates based on copy number thresholds
    """

    def __init__(
        self,
        analyzed_candidates: RecordTypedDF[AnalyzedTnJc2],
        output_dir: Path,
        copy_number_threshold: float = 1.5,
        del_copy_number_threshold: float = 0.3,
        filter_amplicon_length: int = 100,
        force: Optional[bool] = None,
    ):
        self.analyzed_candidates = analyzed_candidates
        self.output_dir = Path(output_dir)
        self.copy_number_threshold = copy_number_threshold
        self.del_copy_number_threshold = del_copy_number_threshold
        self.filter_amplicon_length = filter_amplicon_length
        
        self.isjc2_file = output_dir / "ISJC2.csv"
        self.candidates_file = output_dir / "candidate_amplifications.csv"
        
        super().__init__(
            input_files=[],
            output_files=[self.isjc2_file, self.candidates_file],
            force=force,
        )

    def _calculate_output(self) -> None:
        """Export candidates to CSV."""
        if len(self.analyzed_candidates) == 0:
            info("No candidates to export")
            return
        
        # Convert to DataFrame
        df = self.analyzed_candidates.df.copy()
        
        # Build export DataFrame with correct column names
        export_df = pd.DataFrame()
        
        # Basic columns
        if 'iso_name' in df.columns:
            export_df['isolate'] = df['iso_name']
        if 'ref_name' in df.columns:
            export_df['Reference'] = df['ref_name']
        if 'anc_name' in df.columns:
            export_df['Ancestor'] = df['anc_name']
        
        # Positions
        if 'pos_chr_L' in df.columns and 'pos_chr_R' in df.columns:
            export_df['Positions_in_chromosome'] = (
                df['pos_chr_L'].astype(str) + '-' + df['pos_chr_R'].astype(str)
            )
        
        # Directions
        if 'dir_chr_L' in df.columns and 'dir_chr_R' in df.columns:
            export_df['Direction_in_chromosome'] = (
                df['dir_chr_L'].astype(str) + '/' + df['dir_chr_R'].astype(str)
            )
        
        # Amplicon length
        if 'amplicon_length' in df.columns:
            export_df['amplicon_length'] = df['amplicon_length']
        
        # IS element (format list to string)
        if 'tn_ids' in df.columns:
            export_df['IS_element'] = df['tn_ids'].apply(
                lambda x: ','.join(map(str, x)) if isinstance(x, list) else str(x)
            )
        
        # Coverage columns
        if 'amplicon_coverage' in df.columns:
            export_df['median_copy_number'] = df['amplicon_coverage']
        if 'amplicon_coverage_mode' in df.columns:
            export_df['mode_copy_number'] = df['amplicon_coverage_mode']
        
        # Event columns
        if 'event' in df.columns:
            export_df['event'] = df['event']
        if 'isolate_architecture' in df.columns:
            export_df['isolate_architecture'] = df['isolate_architecture']
        
        # Sort by mode_copy_number descending
        if 'mode_copy_number' in export_df.columns:
            export_df = export_df.sort_values('mode_copy_number', ascending=False)
        
        # Reorder columns
        column_order = [
            'isolate', 'Reference', 'Positions_in_chromosome', 'Direction_in_chromosome',
            'amplicon_length', 'IS_element', 'median_copy_number', 'mode_copy_number',
            'Ancestor', 'event', 'isolate_architecture',
        ]
        available_cols = [c for c in column_order if c in export_df.columns]
        export_df = export_df[available_cols]
        
        # Export ISJC2.csv (all candidates)
        export_df.to_csv(self.isjc2_file, index=False)
        info(f"Exported {len(export_df)} candidates to {self.isjc2_file}")
        
        # Filter candidates for candidate_amplifications.csv
        if 'mode_copy_number' in export_df.columns and 'amplicon_length' in export_df.columns:
            filtered = export_df[
                ((export_df['mode_copy_number'] > self.copy_number_threshold) |
                 (export_df['mode_copy_number'] < self.del_copy_number_threshold)) &
                (export_df['amplicon_length'] > self.filter_amplicon_length)
            ]
            
            filtered.to_csv(self.candidates_file, index=False)
            info(f"Exported {len(filtered)} filtered candidates to {self.candidates_file}")
        else:
            # If columns missing, export empty DataFrame
            export_df.head(0).to_csv(self.candidates_file, index=False)
            info(f"Exported empty candidates file (missing required columns)")

    def _save_output(self, output: None) -> None:
        """Output already saved in _calculate_output."""
        pass

    def load_outputs(self) -> None:
        """Export step has no loadable output."""
        return None
