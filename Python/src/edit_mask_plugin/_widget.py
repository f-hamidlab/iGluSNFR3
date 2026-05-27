from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import tifffile
from qtpy.QtCore import QSettings, Qt
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QSlider,
    QSpinBox,
    QShortcut,
    QVBoxLayout,
    QWidget,
)
from scipy import ndimage as ndi
from scipy.io import loadmat, savemat

MAT_SUFFIX = "_mask_data.mat"


@dataclass
class Record:
    mat_path: Path


class EditMaskWidget(QWidget):
    """Dock widget for browsing folders and editing filtered masks."""

    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        self.launch_root = Path.cwd().resolve()
        self.settings = QSettings("iGluSNFR3", "EditMaskPlugin")
        self.records: List[Record] = []
        self.loaded_files: Dict[Path, Dict[str, object]] = {}
        self.unsaved_files: set[Path] = set()
        self.region_action: Optional[str] = None
        self.scanned_root: Path = self.launch_root
        self.filtered_layer_to_path: Dict[str, Path] = {}
        self.baseline_filtered_masks: Dict[Path, np.ndarray] = {}
        self._syncing_translation = False

        self._build_ui()
        self.viewer.mouse_drag_callbacks.append(self._on_viewer_mouse_drag)
        self.viewer.window.qt_viewer.closeEvent_original = self.viewer.window.qt_viewer.closeEvent
        self.viewer.window.qt_viewer.closeEvent = self._on_window_close

    def _build_ui(self) -> None:
        root_layout = QVBoxLayout(self)

        root_box = QGroupBox("1) Select data directory")
        root_box_layout = QHBoxLayout(root_box)
        self.root_dir_edit = QLineEdit(str(self._initial_root_dir()))
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self._browse_root)
        rescan_button = QPushButton("Scan")
        rescan_button.clicked.connect(self._scan)
        root_box_layout.addWidget(self.root_dir_edit)
        root_box_layout.addWidget(browse_button)
        root_box_layout.addWidget(rescan_button)

        list_box = QGroupBox("2) File list (recursive)")
        list_box_layout = QVBoxLayout(list_box)
        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QAbstractItemView.MultiSelection)
        self.list_widget.itemSelectionChanged.connect(self._on_selection_changed)

        nav_layout = QHBoxLayout()
        add_button = QPushButton("Add to loaded")
        add_button.clicked.connect(self._load_selected_append)
        reset_button = QPushButton("Reset and load")
        reset_button.clicked.connect(self._load_selected_reset)
        select_all_button = QPushButton("Select all")
        select_all_button.clicked.connect(lambda: self._set_list_selection(True))
        deselect_all_button = QPushButton("Deselect all")
        deselect_all_button.clicked.connect(lambda: self._set_list_selection(False))

        unload_button = QPushButton("Unload all")
        unload_button.clicked.connect(self._unload_all)

        nav_layout.addWidget(add_button)
        nav_layout.addWidget(reset_button)
        nav_layout.addWidget(select_all_button)
        nav_layout.addWidget(deselect_all_button)
        nav_layout.addWidget(unload_button)

        list_box_layout.addWidget(self.list_widget)
        list_box_layout.addLayout(nav_layout)

        visibility_box = QGroupBox("3) Layer visibility")
        visibility_layout = QHBoxLayout(visibility_box)
        self.layer_visibility_checks: Dict[str, QCheckBox] = {}

        visibility_defaults = {
            "Projection": True,
            "Filtered image": True,
            "Unfiltered mask": False,
            "Filtered mask": True,
        }

        for suffix, default_visible in visibility_defaults.items():
            checkbox = QCheckBox(suffix)
            checkbox.setChecked(default_visible)
            checkbox.toggled.connect(lambda checked, s=suffix: self._set_layer_type_visibility(s, checked))
            visibility_layout.addWidget(checkbox)
            self.layer_visibility_checks[suffix] = checkbox

        translation_box = QGroupBox("4) Translation")
        translation_layout = QGridLayout(translation_box)

        self.translate_dx_spin = QSpinBox()
        self.translate_dx_spin.setRange(-5000, 5000)
        self.translate_dx_spin.setValue(0)
        self.translate_dy_spin = QSpinBox()
        self.translate_dy_spin.setRange(-5000, 5000)
        self.translate_dy_spin.setValue(0)

        apply_step_button = QPushButton("Apply step")
        apply_step_button.clicked.connect(self._apply_translation_step)
        reset_translation_button = QPushButton("Reset translation")
        reset_translation_button.clicked.connect(self._reset_translation)

        left_button = QPushButton("Left")
        left_button.clicked.connect(lambda: self._nudge_translation(-1, 0))
        right_button = QPushButton("Right")
        right_button.clicked.connect(lambda: self._nudge_translation(1, 0))
        up_button = QPushButton("Up")
        up_button.clicked.connect(lambda: self._nudge_translation(0, -1))
        down_button = QPushButton("Down")
        down_button.clicked.connect(lambda: self._nudge_translation(0, 1))

        self.translation_status = QLabel("Current offset (dx, dy): (0, 0)")

        translation_layout.addWidget(QLabel("Step dx"), 0, 0)
        translation_layout.addWidget(self.translate_dx_spin, 0, 1)
        translation_layout.addWidget(QLabel("Step dy"), 0, 2)
        translation_layout.addWidget(self.translate_dy_spin, 0, 3)
        translation_layout.addWidget(apply_step_button, 0, 4)
        translation_layout.addWidget(reset_translation_button, 0, 5)
        translation_layout.addWidget(left_button, 1, 0)
        translation_layout.addWidget(right_button, 1, 1)
        translation_layout.addWidget(up_button, 1, 2)
        translation_layout.addWidget(down_button, 1, 3)
        shortcut_note = QLabel("Keyboard nudge: Shift+Arrow")
        translation_layout.addWidget(self.translation_status, 2, 0, 1, 6)
        translation_layout.addWidget(shortcut_note, 3, 0, 1, 6)

        edit_box = QGroupBox("5) Edit filtered mask")
        edit_layout = QGridLayout(edit_box)

        self.paint_button = QPushButton("Pencil")
        self.paint_button.setCheckable(True)
        self.paint_button.clicked.connect(lambda checked: self._on_tool_toggled("paint", checked))

        self.erase_button = QPushButton("Eraser")
        self.erase_button.setCheckable(True)
        self.erase_button.clicked.connect(lambda checked: self._on_tool_toggled("erase", checked))

        self.add_region_button = QPushButton("Add region")
        self.add_region_button.setCheckable(True)
        self.add_region_button.clicked.connect(lambda checked: self._on_tool_toggled("add_region", checked))

        self.delete_region_button = QPushButton("Delete region")
        self.delete_region_button.setCheckable(True)
        self.delete_region_button.clicked.connect(lambda checked: self._on_tool_toggled("delete_region", checked))

        brush_label = QLabel("Brush size")
        self.brush_slider = QSlider(Qt.Horizontal)
        self.brush_slider.setMinimum(1)
        self.brush_slider.setMaximum(100)
        self.brush_slider.setValue(10)
        self.brush_slider.valueChanged.connect(self._on_brush_size_changed)
        self.brush_size_value = QSpinBox()
        self.brush_size_value.setMinimum(1)
        self.brush_size_value.setMaximum(100)
        self.brush_size_value.setValue(self.brush_slider.value())
        self.brush_size_value.valueChanged.connect(self._on_brush_size_spin_changed)

        undo_note = QLabel("Tip: Undo edits with Ctrl+Z")

        edit_layout.addWidget(self.paint_button, 0, 0)
        edit_layout.addWidget(self.erase_button, 0, 1)
        edit_layout.addWidget(self.add_region_button, 0, 2)
        edit_layout.addWidget(self.delete_region_button, 0, 3)
        edit_layout.addWidget(brush_label, 1, 0)
        edit_layout.addWidget(self.brush_slider, 1, 1, 1, 2)
        edit_layout.addWidget(self.brush_size_value, 1, 3)
        edit_layout.addWidget(undo_note, 2, 0, 1, 4)

        save_box = QGroupBox("6) Save")
        save_layout = QVBoxLayout(save_box)
        save_note = QLabel("Save overwrites *_mask_data.mat and writes *_binary.tif")
        self.unsaved_label = QLabel("No unsaved changes")
        self.unsaved_label.setStyleSheet("color: green;")

        save_btn_layout = QHBoxLayout()
        save_active_button = QPushButton("Save active")
        save_active_button.clicked.connect(self._save_active)
        save_all_button = QPushButton("Save all modified")
        save_all_button.clicked.connect(self._save_all_modified)
        save_btn_layout.addWidget(save_active_button)
        save_btn_layout.addWidget(save_all_button)

        save_layout.addWidget(save_note)
        save_layout.addWidget(self.unsaved_label)
        save_layout.addLayout(save_btn_layout)

        self.status_label = QLabel("Ready")

        root_layout.addWidget(root_box)
        root_layout.addWidget(list_box)
        root_layout.addWidget(visibility_box)
        root_layout.addWidget(translation_box)
        root_layout.addWidget(edit_box)
        root_layout.addWidget(save_box)
        root_layout.addWidget(self.status_label)

        self._setup_translation_shortcuts()
        self._scan()

    def _setup_translation_shortcuts(self) -> None:
        shortcut_bindings = [
            ("Shift+Left", -1, 0),
            ("Shift+Right", 1, 0),
            ("Shift+Up", 0, -1),
            ("Shift+Down", 0, 1),
        ]

        self._translation_shortcuts: List[QShortcut] = []
        for keyseq, dx, dy in shortcut_bindings:
            shortcut = QShortcut(keyseq, self.viewer.window.qt_viewer)
            shortcut.setContext(Qt.ApplicationShortcut)
            shortcut.activated.connect(lambda dx=dx, dy=dy: self._nudge_translation(dx, dy))
            self._translation_shortcuts.append(shortcut)

    def _initial_root_dir(self) -> Path:
        saved_root = self.settings.value("last_root_dir", "", type=str)
        if saved_root:
            saved_path = Path(saved_root).expanduser()
            if saved_path.exists():
                return saved_path.resolve()
        return self.launch_root

    def _store_last_root(self, root: Path) -> None:
        self.settings.setValue("last_root_dir", str(root))

    def _browse_root(self) -> None:
        selected = QFileDialog.getExistingDirectory(self, "Select root folder", self.root_dir_edit.text())
        if selected:
            self.root_dir_edit.setText(selected)
            self._scan()

    def _scan(self) -> None:
        root = Path(self.root_dir_edit.text()).expanduser().resolve()
        if not root.exists():
            self.root_dir_edit.setText(str(self.launch_root))
            self.status_label.setText("Root folder does not exist; reverted to Napari launch directory")
            return

        self.scanned_root = root
        self._store_last_root(root)
        self.records = self._discover_records(root)
        self.list_widget.clear()

        for rec in self.records:
            self.list_widget.addItem(str(rec.mat_path.relative_to(root)))

        if self.records:
            self.status_label.setText(f"Found {len(self.records)} mask-data files")
        else:
            self.status_label.setText("No *_mask_data.mat files found")

    def _discover_records(self, root: Path) -> List[Record]:
        records = [Record(mat_path=p) for p in root.rglob(f"*{MAT_SUFFIX}")]
        records.sort(key=lambda r: str(r.mat_path))
        return records

    def _dataset_key(self, mat_path: Path) -> str:
        try:
            rel = mat_path.relative_to(self.scanned_root)
            label = str(rel)
        except ValueError:
            label = str(mat_path)

        if label.endswith(MAT_SUFFIX):
            label = label[: -len(MAT_SUFFIX)]
        elif label.endswith(".mat"):
            label = label[:-4]

        return label

    def _layer_name(self, mat_path: Path, suffix: str) -> str:
        return f"{self._dataset_key(mat_path)} | {suffix}"

    def _display_path(self, mat_path: Path) -> str:
        try:
            return str(mat_path.relative_to(self.scanned_root))
        except ValueError:
            return str(mat_path)

    def _translation_setting_key(self, mat_path: Path) -> str:
        return str(mat_path.resolve()).replace("/", "__")

    def _save_translation_offset(self, mat_path: Path, dx: float, dy: float) -> None:
        key = self._translation_setting_key(mat_path)
        self.settings.setValue(f"translation_offsets/{key}/dx", float(dx))
        self.settings.setValue(f"translation_offsets/{key}/dy", float(dy))

    def _load_translation_offset(self, mat_path: Path) -> tuple[float, float]:
        key = self._translation_setting_key(mat_path)
        dx = self.settings.value(f"translation_offsets/{key}/dx", 0.0, type=float)
        dy = self.settings.value(f"translation_offsets/{key}/dy", 0.0, type=float)
        return float(dx), float(dy)

    def _is_layer_type_visible(self, suffix: str) -> bool:
        checkbox = self.layer_visibility_checks.get(suffix)
        return checkbox.isChecked() if checkbox is not None else True

    def _set_layer_type_visibility(self, suffix: str, visible: bool) -> None:
        for mat_path in self.loaded_files:
            layer_name = self._layer_name(mat_path, suffix)
            if layer_name in self.viewer.layers:
                self.viewer.layers[layer_name].visible = visible

    def _on_selection_changed(self) -> None:
        self._update_translation_status()

    def _set_list_selection(self, select: bool) -> None:
        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)
            item.setSelected(select)

    def _load_selected_append(self) -> None:
        self._load_selected(append=True)

    def _load_selected_reset(self) -> None:
        if self.unsaved_files:
            unsaved_list = ", ".join(self._display_path(p) for p in sorted(self.unsaved_files))
            reply = self._confirm_dialog(
                f"Unsaved changes in: {unsaved_list}\n\nDiscard and reset?",
                "Reset and Load",
            )
            if not reply:
                return
        self._load_selected(append=False)

    def _confirm_dialog(self, message: str, title: str = "Confirm") -> bool:
        msgbox = QMessageBox(self)
        msgbox.setWindowTitle(title)
        msgbox.setText(message)
        msgbox.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)
        return msgbox.exec() == QMessageBox.Yes

    def _load_selected(self, append: bool = True) -> None:
        selected_indices = sorted({i.row() for i in self.list_widget.selectedIndexes()})
        if not selected_indices:
            self.status_label.setText("Select one or more files")
            return

        if not append:
            if not self._unload_all(force=True):
                return

        loaded_count = 0
        for idx in selected_indices:
            if idx < 0 or idx >= len(self.records):
                continue
            rec = self.records[idx]
            if rec.mat_path in self.loaded_files:
                continue

            mat = loadmat(rec.mat_path)
            mat_data = {k: v for k, v in mat.items() if not k.startswith("__")}

            required = ["im_denoised", "im_filtered", "binary_mask", "filtered_mask"]
            missing = [k for k in required if k not in mat_data]
            if missing:
                self.status_label.setText(f"Skipping {rec.mat_path.name}: missing {', '.join(missing)}")
                continue

            self.loaded_files[rec.mat_path] = mat_data
            self._create_layers_for_file(rec.mat_path, mat_data)
            loaded_count += 1

        self._update_unsaved_label()
        self.status_label.setText(f"Loaded {loaded_count} file(s)")

    def _create_layers_for_file(self, mat_path: Path, mat_data: Dict[str, object]) -> None:
        projection = self._to_image(mat_data["im_denoised"])
        filtered_image = self._to_image(mat_data["im_filtered"])
        unfiltered_mask = self._to_mask(mat_data["binary_mask"])
        filtered_mask = self._to_mask(mat_data["filtered_mask"])
        self.baseline_filtered_masks[mat_path] = (filtered_mask > 0).astype(np.uint8)
        unfiltered_mask = (unfiltered_mask > 0).astype(np.uint8) * 2

        self._set_or_replace_image(self._layer_name(mat_path, "Projection"), projection, colormap="gray", opacity=0.8)
        self._set_or_replace_image(
            self._layer_name(mat_path, "Filtered image"), filtered_image, colormap="green", opacity=0.8
        )
        self._set_or_replace_labels(
            self._layer_name(mat_path, "Unfiltered mask"), unfiltered_mask, editable=False, opacity=0.35
        )
        if self._layer_name(mat_path, "Unfiltered mask") in self.viewer.layers:
            self.viewer.layers[self._layer_name(mat_path, "Unfiltered mask")].visible = False

        self._set_or_replace_labels(
            self._layer_name(mat_path, "Filtered mask"), filtered_mask, editable=True, opacity=0.45
        )
        for suffix in ["Projection", "Filtered image", "Unfiltered mask", "Filtered mask"]:
            layer_name = self._layer_name(mat_path, suffix)
            if layer_name in self.viewer.layers:
                self.viewer.layers[layer_name].visible = self._is_layer_type_visible(suffix)

        self.filtered_layer_to_path[self._layer_name(mat_path, "Filtered mask")] = mat_path
        filtered_layer = self.viewer.layers[self._layer_name(mat_path, "Filtered mask")]
        self._connect_unsaved_tracking(filtered_layer, mat_path)
        projection_layer = self.viewer.layers[self._layer_name(mat_path, "Projection")]
        projection_layer.events.translate.connect(lambda _e, p=mat_path: self._sync_translation_from_projection(p))

        self._set_brush_size(self.brush_slider.value())
        self._activate_tool("paint")
        saved_dx, saved_dy = self._load_translation_offset(mat_path)
        if saved_dx != 0.0 or saved_dy != 0.0:
            self._set_translation_offset(mat_path, saved_dx, saved_dy)
        self._update_translation_status()

    def _set_or_replace_image(self, name: str, arr: np.ndarray, colormap: str, opacity: float) -> None:
        if name in self.viewer.layers:
            layer = self.viewer.layers[name]
            layer.data = arr
            layer.colormap = colormap
            layer.opacity = opacity
        else:
            self.viewer.add_image(arr, name=name, colormap=colormap, opacity=opacity, blending="additive")

    def _set_or_replace_labels(self, name: str, labels_data: np.ndarray, editable: bool, opacity: float) -> None:
        if name in self.viewer.layers:
            layer = self.viewer.layers[name]
            layer.data = labels_data
            layer.opacity = opacity
        else:
            layer = self.viewer.add_labels(labels_data, name=name, opacity=opacity)

        layer.editable = editable
        if editable:
            layer.selected_label = 1
            layer.mode = "paint"

    def _connect_unsaved_tracking(self, layer, mat_path: Path) -> None:
        # Different napari edit paths may emit different events (in-place paint/erase vs full data replace).
        for event_name in ["data", "set_data", "labels_update", "paint", "erase"]:
            event_emitter = getattr(layer.events, event_name, None)
            if event_emitter is not None:
                event_emitter.connect(lambda _e, p=mat_path: self._mark_unsaved(p))

    def _editable_layer(self):
        active = self.viewer.layers.selection.active
        if active is not None and active.name in self.filtered_layer_to_path:
            return active

        for layer_name in self.filtered_layer_to_path:
            if layer_name in self.viewer.layers:
                return self.viewer.layers[layer_name]

        self.status_label.setText("Load a record with filtered mask layers first")
        return None

    def _on_tool_toggled(self, tool: str, checked: bool) -> None:
        if checked:
            self._activate_tool(tool)
        else:
            self._activate_tool("paint")

    def _activate_tool(self, tool: str) -> None:
        layer = self._editable_layer()
        if layer is None:
            return

        buttons = {
            "paint": self.paint_button,
            "erase": self.erase_button,
            "add_region": self.add_region_button,
            "delete_region": self.delete_region_button,
        }

        for key, button in buttons.items():
            button.blockSignals(True)
            button.setChecked(key == tool)
            button.blockSignals(False)

        if tool == "paint":
            self.region_action = None
            layer.mode = "paint"
            layer.selected_label = 1
            self.status_label.setText("Pencil tool enabled")
        elif tool == "erase":
            self.region_action = None
            layer.mode = "erase"
            self.status_label.setText("Eraser tool enabled")
        elif tool == "add_region":
            self.region_action = "add"
            layer.mode = "pan_zoom"
            self.status_label.setText("Add region enabled: click a connected region in Unfiltered mask")
        elif tool == "delete_region":
            self.region_action = "delete"
            layer.mode = "pan_zoom"
            self.status_label.setText("Delete region enabled: click a connected region in Unfiltered mask")

    def _mark_unsaved(self, mat_path: Path) -> None:
        if mat_path not in self.loaded_files:
            return

        layer_name = self._layer_name(mat_path, "Filtered mask")
        if layer_name not in self.viewer.layers:
            return

        current_mask = (np.asarray(self.viewer.layers[layer_name].data) > 0).astype(np.uint8)
        baseline = self.baseline_filtered_masks.get(mat_path)

        if baseline is None:
            self.baseline_filtered_masks[mat_path] = current_mask.copy()
            self.unsaved_files.discard(mat_path)
        elif current_mask.shape == baseline.shape and np.array_equal(current_mask, baseline):
            self.unsaved_files.discard(mat_path)
        else:
            self.unsaved_files.add(mat_path)

        self._update_unsaved_label()

    def _update_unsaved_label(self) -> None:
        if not self.unsaved_files:
            self.unsaved_label.setText("No unsaved changes")
            self.unsaved_label.setStyleSheet("color: green;")
        else:
            unsaved_names = ", ".join(self._display_path(p) for p in sorted(self.unsaved_files))
            self.unsaved_label.setText(f"Modified but not saved: {unsaved_names}")
            self.unsaved_label.setStyleSheet("color: red;")

    def _on_brush_size_changed(self, size: int) -> None:
        self.brush_size_value.blockSignals(True)
        self.brush_size_value.setValue(size)
        self.brush_size_value.blockSignals(False)
        self._set_brush_size(size)

    def _on_brush_size_spin_changed(self, size: int) -> None:
        self.brush_slider.blockSignals(True)
        self.brush_slider.setValue(size)
        self.brush_slider.blockSignals(False)
        self._set_brush_size(size)

    def _set_brush_size(self, size: int) -> None:
        layer = self._editable_layer()
        if layer is None:
            return
        layer.brush_size = float(size)

    def _get_active_mat_path(self) -> Optional[Path]:
        active = self.viewer.layers.selection.active
        if active is not None:
            return self.filtered_layer_to_path.get(active.name)
        return None

    def _get_active_dataset_path(self) -> Optional[Path]:
        active = self.viewer.layers.selection.active
        if active is None:
            return None

        # Prefer direct filtered-mask mapping when available.
        mapped = self.filtered_layer_to_path.get(active.name)
        if mapped is not None:
            return mapped

        # Fall back to matching any known layer name for loaded datasets.
        for mat_path in self.loaded_files:
            for suffix in ["Projection", "Filtered image", "Unfiltered mask", "Filtered mask"]:
                if active.name == self._layer_name(mat_path, suffix):
                    return mat_path
        return None

    def _dataset_layers(self, mat_path: Path) -> Dict[str, object]:
        out: Dict[str, object] = {}
        for suffix in ["Projection", "Filtered image", "Unfiltered mask", "Filtered mask"]:
            name = self._layer_name(mat_path, suffix)
            if name in self.viewer.layers:
                out[suffix] = self.viewer.layers[name]
        return out

    def _set_translation_offset(self, mat_path: Path, dx: float, dy: float) -> None:
        layers = self._dataset_layers(mat_path)
        for layer in layers.values():
            tr = list(layer.translate)
            if len(tr) >= 2:
                tr[-2] = float(dy)
                tr[-1] = float(dx)
                layer.translate = tuple(tr)
        self._save_translation_offset(mat_path, float(dx), float(dy))
        self._update_translation_status()

    def _apply_translation_delta(self, mat_path: Path, dx: float, dy: float) -> None:
        current_dx, current_dy = self._translation_offset(mat_path)
        self._set_translation_offset(mat_path, current_dx + float(dx), current_dy + float(dy))

    def _translation_offset(self, mat_path: Path) -> tuple[float, float]:
        name = self._layer_name(mat_path, "Projection")
        if name not in self.viewer.layers:
            return (0.0, 0.0)
        tr = self.viewer.layers[name].translate
        if len(tr) < 2:
            return (0.0, 0.0)
        return (float(tr[-1]), float(tr[-2]))

    def _update_translation_status(self) -> None:
        mat_path = self._get_active_dataset_path()
        if mat_path is None:
            self.translation_status.setText("Current offset (dx, dy): (0, 0)")
            return
        dx, dy = self._translation_offset(mat_path)
        self.translation_status.setText(f"Current offset (dx, dy): ({dx:.1f}, {dy:.1f})")

    def _apply_translation_step(self) -> None:
        mat_path = self._get_active_dataset_path()
        if mat_path is None:
            self.status_label.setText("Select any layer from a loaded dataset first")
            return
        self._apply_translation_delta(mat_path, self.translate_dx_spin.value(), self.translate_dy_spin.value())

    def _nudge_translation(self, dx: int, dy: int) -> None:
        mat_path = self._get_active_dataset_path()
        if mat_path is None:
            self.status_label.setText("Select any layer from a loaded dataset first")
            return
        self._apply_translation_delta(mat_path, dx, dy)

    def _reset_translation(self) -> None:
        mat_path = self._get_active_dataset_path()
        if mat_path is None:
            self.status_label.setText("Select any layer from a loaded dataset first")
            return
        self._set_translation_offset(mat_path, 0.0, 0.0)

    def _sync_translation_from_projection(self, mat_path: Path) -> None:
        if self._syncing_translation:
            return
        name = self._layer_name(mat_path, "Projection")
        if name not in self.viewer.layers:
            return

        self._syncing_translation = True
        try:
            tr = tuple(self.viewer.layers[name].translate)
            for suffix in ["Filtered image", "Unfiltered mask", "Filtered mask"]:
                lname = self._layer_name(mat_path, suffix)
                if lname in self.viewer.layers:
                    self.viewer.layers[lname].translate = tr
            if len(tr) >= 2:
                self._save_translation_offset(mat_path, float(tr[-1]), float(tr[-2]))
        finally:
            self._syncing_translation = False
        self._update_translation_status()

    def _on_viewer_mouse_drag(self, viewer, event) -> None:
        if event.type != "mouse_press":
            return
        if self.region_action not in {"add", "delete"}:
            return
        self._apply_region_action_from_unfiltered()

    def _apply_region_action_from_unfiltered(self) -> None:
        layer = self._editable_layer()
        if layer is None:
            return

        mat_path = self._get_active_mat_path()
        if mat_path is None:
            self.status_label.setText("Select a filtered mask layer first")
            return

        source_name = self._layer_name(mat_path, "Unfiltered mask")
        if source_name not in self.viewer.layers:
            self.status_label.setText("Unfiltered mask layer is missing")
            return

        source_layer = self.viewer.layers[source_name]
        source_data = np.asarray(source_layer.data)
        if source_data.ndim != 2:
            self.status_label.setText("Region actions currently support 2D masks")
            return

        try:
            world_pos = self.viewer.cursor.position
            source_coord = tuple(int(round(c)) for c in source_layer.world_to_data(world_pos))
        except Exception:
            self.status_label.setText("Move mouse over a region in Unfiltered mask")
            return

        if len(source_coord) < 2:
            self.status_label.setText("Invalid cursor position")
            return

        sr, sc = source_coord[-2], source_coord[-1]
        if sr < 0 or sc < 0 or sr >= source_data.shape[0] or sc >= source_data.shape[1] or source_data[sr, sc] == 0:
            self.status_label.setText("Cursor is not on a region in Unfiltered mask")
            return

        source_components, _ = ndi.label(source_data > 0)
        component_id = source_components[sr, sc]
        if component_id == 0:
            self.status_label.setText("No source component under cursor")
            return

        component_mask = source_components == component_id

        data = np.asarray(layer.data).copy()
        if data.ndim != 2:
            self.status_label.setText("Region actions currently support 2D masks")
            return

        if self.region_action == "add":
            data[component_mask] = 1
            action_text = "Added"
        else:
            data[component_mask] = 0
            action_text = "Deleted"

        layer.data = data.astype(np.uint8)
        self._mark_unsaved(mat_path)
        self.status_label.setText(f"{action_text} connected region from unfiltered mask")

    def _save_active(self) -> None:
        mat_path = self._get_active_mat_path()
        if mat_path is None:
            self.status_label.setText("No active filtered mask layer")
            return
        self._save_single_file(mat_path)
        self.unsaved_files.discard(mat_path)
        self._update_unsaved_label()

    def _save_all_modified(self) -> None:
        if not self.unsaved_files:
            self.status_label.setText("No unsaved changes")
            return

        for mat_path in list(self.unsaved_files):
            self._save_single_file(mat_path)
        self.unsaved_files.clear()
        self._update_unsaved_label()

    def _save_single_file(self, mat_path: Path) -> None:
        if mat_path not in self.loaded_files:
            self.status_label.setText(f"File not loaded: {mat_path.name}")
            return

        layer_name = self._layer_name(mat_path, "Filtered mask")
        if layer_name not in self.viewer.layers:
            self.status_label.setText(f"Layer not found: {layer_name}")
            return

        layer = self.viewer.layers[layer_name]
        data = (np.asarray(layer.data) > 0).astype(np.uint8)

        mat_out = dict(self.loaded_files[mat_path])
        mat_out["filtered_mask"] = data.astype(bool)

        savemat(mat_path, mat_out)
        stem = mat_path.name
        if stem.endswith(MAT_SUFFIX):
            stem = stem[: -len(MAT_SUFFIX)]
        binary_path = mat_path.with_name(f"{stem}_binary.tif")
        tifffile.imwrite(binary_path, data * 255)

        self.baseline_filtered_masks[mat_path] = data.copy()

        self.status_label.setText(f"Saved: {mat_path.name}")

    def _unload_all(self, force: bool = False) -> bool:
        if self.unsaved_files and not force:
            unsaved_list = ", ".join(self._display_path(p) for p in sorted(self.unsaved_files))
            reply = self._confirm_dialog(
                f"Unsaved changes in: {unsaved_list}\n\nDiscard and unload?",
                "Unload All",
            )
            if not reply:
                return False

        for mat_path in list(self.loaded_files.keys()):
            for suffix in ["Projection", "Filtered image", "Unfiltered mask", "Filtered mask"]:
                layer_name = self._layer_name(mat_path, suffix)
                if layer_name in self.viewer.layers:
                    self.viewer.layers.remove(layer_name)

        self.loaded_files.clear()
        self.unsaved_files.clear()
        self.baseline_filtered_masks.clear()
        self.filtered_layer_to_path.clear()
        self._update_unsaved_label()
        self.status_label.setText("All files unloaded")
        return True

    def _on_window_close(self, event) -> None:
        if self.unsaved_files:
            unsaved_list = ", ".join(self._display_path(p) for p in sorted(self.unsaved_files))
            reply = self._confirm_dialog(
                f"Unsaved changes in: {unsaved_list}\n\nDiscard and close?",
                "Close Window",
            )
            if not reply:
                event.ignore()
                return
        self.viewer.window.qt_viewer.closeEvent_original(event)

    def _to_image(self, arr: object) -> np.ndarray:
        img = np.asarray(arr)
        img = np.squeeze(img)

        if img.ndim > 2:
            while img.ndim > 2:
                img = img.max(axis=0)
        return img

    def _to_mask(self, arr: object) -> np.ndarray:
        mask = np.asarray(arr)
        mask = np.squeeze(mask)
        if mask.ndim > 2:
            while mask.ndim > 2:
                mask = mask.max(axis=0)
        return (mask > 0).astype(np.uint8)
