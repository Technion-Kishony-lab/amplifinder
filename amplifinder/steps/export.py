"""Step 14: Export results to YAML."""

from pathlib import Path
from typing import Optional, Dict, Any
from dataclasses import asdict
import yaml

from amplifinder.data_types import RecordTypedDf, ClassifiedTnJc2, Genome
from amplifinder.steps.base import OutputStep
from amplifinder.steps.read_length import ReadLengths


class ExportTnJc2Step(OutputStep[Dict[str, Any]]):
    """Export classified amplicons to YAML."""

    def __init__(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        genome: Genome,
        output_dir: Path,
        ref_name: str,
        iso_name: str,
        read_lengths: ReadLengths,
        anc_name: Optional[str] = None,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self.genome = genome
        self.ref_name = ref_name
        self.iso_name = iso_name
        self.anc_name = anc_name
        self.read_lengths = read_lengths

        self.yaml_file = output_dir / "amplifications.yaml"

        super().__init__(output_files=[self.yaml_file], force=force)

    def _calculate_output(self) -> Dict[str, Any]:
        """Build export structure."""
        amplicons = []
        for tnjc2 in self.classified_tnjc2s:
            amplicons.append({
                'positions': f"{tnjc2.left}-{tnjc2.right}",
                'span_origin': tnjc2.span_origin,
                'length': tnjc2.amplicon_length,
                'IS_elements': tnjc2.tn_ids,
                'chosen_tn_id': tnjc2.chosen_tn_id,
                'copy_number': tnjc2.copy_number,
                'architecture': tnjc2.iso_architecture.description,
                'descriptors': [event_descriptor.value for event_descriptor in tnjc2.event_descriptors],
            })

        return {
            'sample': {
                'isolate': self.iso_name,
                'ancestor': self.anc_name,
                'reference': self.ref_name,
                'read_lengths': asdict(self.read_lengths),
            },
            'amplicons': amplicons,
        }

    def _save_output(self, output: Dict[str, Any]) -> None:
        """Write YAML file."""
        with open(self.yaml_file, 'w') as f:
            yaml.dump(output, f, default_flow_style=False, sort_keys=False)

    def report_output_message(self, output: Dict[str, Any]) -> Optional[str]:
        return f"Exported {len(output['amplicons'])} amplicons to {self.yaml_file}"
