from PySide6.QtCore import Qt, QAbstractListModel, QSize, QRect
from PySide6.QtGui import QPainter, QPixmap, QFont, QColor
from PySide6.QtWidgets import (QMdiSubWindow, QWidget, QGraphicsDropShadowEffect, QStyleOption, QStyle, QVBoxLayout, QListView, QFrame,
                                QStyledItemDelegate, QLabel)


class Delegate(QStyledItemDelegate):
    def __init__(self, height=None):
        super(Delegate, self).__init__()
        if height is None:
            self._height = 45
        else:
            self._height = height

    def paint(self, painter, option, index):
        super(Delegate, self).paint(painter, option, index)

        # HOVER
        if option.state & QStyle.State_MouseOver:
            painter.fillRect(option.rect, QColor("#F1F1F1"))
        else:
            painter.fillRect(option.rect, Qt.transparent)

        # SELECTED
        if option.state & QStyle.State_Selected:
            painter.fillRect(option.rect, QColor("#F1F1F1"))

        # DRAW ICON
        icon = QPixmap()
        icon.load(index.data()[1])
        icon = icon.scaled(24, 24, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)

        left = 10 # margin left
        icon_pos = QRect(left, ((self._height - icon.height()) / 2) + option.rect.y(),
                icon.width(), icon.height())
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.SmoothPixmapTransform)
        painter.drawPixmap(icon_pos, icon)

        # DRAW TEXT
        font = QFont("Roboto Black", 12)
        text_pos = QRect((left * 2) + icon.width(), option.rect.y(),
                option.rect.width(), option.rect.height())
        painter.setFont(font)
        painter.setPen(Qt.black)
        painter.drawText(text_pos, Qt.AlignVCenter, index.data()[0])



    def sizeHint(self, option, index):
        return QSize(0, self._height)


class Model(QAbstractListModel):
    def __init__(self, data):
        super(Model, self).__init__()
        self._data = data

    def rowCount(self, index):
        return len(self._data)

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid() and role == Qt.DisplayRole:
            return self._data[index.row()]


class ListView(QListView):
    def __init__(self):
        super(ListView, self).__init__()
        self.setMouseTracking(True)
        
    def mouseMoveEvent(self, event):
        # CHANGE CURSOR HOVERING
        if self.indexAt(event.pos()).row() >= 0:
            self.setCursor(Qt.PointingHandCursor)
        else:
            self.setCursor(Qt.ArrowCursor)


class SideMenuWidget(QWidget):
    def __init__(self, data):
        super().__init__()
        self.layout = QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        self.listview = ListView()
        self.listview.setModel(Model(data))
        self.listview.setItemDelegate(Delegate())
        self.layout.addWidget(self.listview)

        self.setStyleSheet("background: white;")

    def sizeHint(self):
        size_hint = self.listview.sizeHint()
        size_hint.setWidth(200)
        return size_hint
