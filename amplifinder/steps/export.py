"""Step 14: Export results to YAML."""

from pathlib import Path
from typing import Optional, Dict, Any
from dataclasses import asdict
import yaml
import json

from amplifinder.data_types import RecordTypedDf, ClassifiedTnJc2, CoveredTnJc2, RefTn, BaseEvent
from amplifinder.steps.base import OutputStep
from amplifinder.steps.read_length import ReadLengths
from amplifinder.utils.json_utils import compact_short_lists


class ExportTnJc2Step(OutputStep[Dict[str, Any]]):
    """Export classified amplicons to YAML."""

    def __init__(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        linked_tnjc2s: RecordTypedDf[CoveredTnJc2],
        output_dir: Path,
        ref_name: str,
        iso_name: str,
        read_lengths: ReadLengths,
        ref_tns: RecordTypedDf[RefTn],
        anc_name: Optional[str] = None,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self.linked_tnjc2s = linked_tnjc2s
        self.ref_name = ref_name
        self.iso_name = iso_name
        self.anc_name = anc_name
        self.read_lengths = read_lengths
        self.ref_tns = ref_tns

        self.yaml_file = output_dir / "summary.yaml"
        self.json_file = output_dir / "summary.json"

        super().__init__(output_files=[self.yaml_file, self.json_file], force=force)

    def _tn_ids_to_names(self, tnjc2: CoveredTnJc2 | ClassifiedTnJc2) -> tuple[list[str], Optional[str]]:
        """Map TN IDs to names."""
        tn_names = [self.ref_tns.df.loc[tn_id, 'tn_name'] for tn_id in tnjc2.tn_ids]
        chosen_tn_name = self.ref_tns.df.loc[tnjc2.chosen_tn_id, 'tn_name'] if tnjc2.chosen_tn_id else None
        return tn_names, chosen_tn_name

    def _calculate_output(self) -> Dict[str, Any]:
        """Build export structure."""
        amplicons = []

        # Export classified amplicons (went through full pipeline)
        for tnjc2 in self.classified_tnjc2s:
            tn_names, chosen_tn_name = self._tn_ids_to_names(tnjc2)

            amplicons.append({
                'positions': f"{tnjc2.left}-{tnjc2.right}",
                'span_origin': tnjc2.span_origin,
                'length': tnjc2.amplicon_length,
                'possible_ISs': tn_names,
                'chosen_IS': chosen_tn_name,
                'copy_number': round(tnjc2.copy_number, 1),
                'architecture': tnjc2.iso_architecture.description,
                'ancestor_architecture': tnjc2.anc_architecture.description if tnjc2.anc_architecture else None,
                'descriptors': [event_descriptor.value for event_descriptor in tnjc2.event_descriptors],
                'left_is_ref_tn': tnjc2.tnjc_left.is_ref_tn_junction(),
                'right_is_ref_tn': tnjc2.tnjc_right.is_ref_tn_junction(),
            })

        # Collect TRANSPOSITION events (filtered out from main pipeline)
        transpositions = []
        for tnjc2 in self.linked_tnjc2s:
            if tnjc2.base_event != BaseEvent.TRANSPOSITION:
                continue
            tn_names, chosen_tn_name = self._tn_ids_to_names(tnjc2)
            transpositions.append({
                'positions': f"{tnjc2.left}-{tnjc2.right}",
                # Not using amplicon_length because if neg we get genome len
                'length': tnjc2.right - tnjc2.left,
                'possible_ISs': tn_names,
            })

        return {
            'sample': {
                'isolate': self.iso_name,
                'ancestor': self.anc_name,
                'reference': self.ref_name,
                'read_lengths': asdict(self.read_lengths),
            },
            'amplicons': amplicons,
            'transpositions': transpositions,
        }

    def _save_output(self, output: Dict[str, Any]) -> None:
        """Write YAML and JSON files."""
        with open(self.yaml_file, 'w') as f:
            yaml.dump(output, f, default_flow_style=False, sort_keys=False)
        json_str = json.dumps(output, indent=2, sort_keys=False)
        json_str = compact_short_lists(json_str)
        with open(self.json_file, 'w') as f:
            f.write(json_str)

    def report_output_message(self, output: Dict[str, Any]) -> Optional[str]:
        return (f"Exported {len(output['amplicons'])} amplicons, "
                f"{len(output['transpositions'])} transpositions")