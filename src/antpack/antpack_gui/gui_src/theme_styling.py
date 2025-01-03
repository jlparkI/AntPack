"""Basic tool for applying a stylesheet to a PySide app.
The code is from qt_material but is reproduced here (with
some slight simplifications) so that we do not need to
maintain that package as a dependency. This is currently
in development and will be used in a future release."""
import os
import shutil
import sys
import platform
from xml.dom.minidom import parse
import jinja2
from PySide6.QtGui import QPalette, QGuiApplication



def use_css_to_style(
    app,
    theme='',
    parent='theme',
    css_file=None,
):
    """"""
    stylesheet = build_stylesheet(theme, parent)
    if stylesheet is None:
        return

    if css_file and os.path.exists(css_file):
        with open(css_file) as file:
            stylesheet += file.read().format(**os.environ)

    try:
        app.setStyleSheet(stylesheet)
    except:
        app.style_sheet = stylesheet



def build_stylesheet(
    theme='',
    parent='theme',
    template=None,
    export=False,
):
    """"""
    if not export:
        try:
            add_fonts()
        except Exception as e:
            logging.warning(e)

    theme = get_theme_colors()
    if theme is None:
        return None

    set_icons_theme(theme, parent=parent)

    # Render custom template
    if os.path.exists(template):
        parent, template = os.path.split(template)
        loader = jinja2.FileSystemLoader(parent)
        env = jinja2.Environment(autoescape=False, loader=loader)
        env.filters['opacity'] = opacity
        env.filters['density'] = density
        stylesheet = env.get_template(template)
    else:
        env = jinja2.Environment(autoescape=False, loader=jinja2.BaseLoader)
        env.filters['opacity'] = opacity
        env.filters['density'] = density
        stylesheet.from_string(template)

    theme.setdefault('icon', None)
    theme.setdefault('font_family', 'Roboto')
    theme.setdefault('danger', '#dc3545')
    theme.setdefault('warning', '#ffc107')
    theme.setdefault('success', '#17a2b8')
    theme.setdefault('density_scale', '0')
    theme.setdefault('button_shape', 'default')

    theme.update(extra)

    default_palette = QGuiApplication.palette()
    color = QColor(
        *[
            int(theme['primaryColor'][i : i + 2], 16)
            for i in range(1, 6, 2)
        ]
        + [92]
    )

    default_palette.set_color(QPalette.ColorRole.Text, color)
    QGuiApplication.set_palette(default_palette)

    environ = {
        'linux': platform.system() == 'Linux',
        'windows': platform.system() == 'Windows',
        'darwin': platform.system() == 'Darwin',
        'pyqt5': 'PyQt5' in sys.modules,
        'pyqt6': 'PyQt6' in sys.modules,
        'pyside2': 'PySide2' in sys.modules,
        'pyside6': 'PySide6' in sys.modules,
    }

    environ.update(theme)
    return stylesheet.render(environ)



def get_theme_colors():
    theme = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                "theme_template.xml")

    document = parse(theme)
    theme = {
        child.getAttribute('name'): child.firstChild.nodeValue
        for child in document.getElementsByTagName('color')
    }

    for k in theme:
        os.environ[str(k)] = theme[k]

    for color in [
        'primaryColor',
        'primaryLightColor',
        'secondaryColor',
        'secondaryLightColor',
        'secondaryDarkColor',
        'primaryTextColor',
        'secondaryTextColor',
    ]:
        os.environ[f'QTMATERIAL_{color.upper()}'] = theme[color]
    os.environ['QTMATERIAL_THEME'] = "Default"

    return theme



def set_icons_theme(theme, parent='theme'):
    """"""
    source = os.path.join(os.path.dirname(__file__), 'resources', 'source')
    resources = ResourceGenerator(
        primary=theme['primaryColor'],
        secondary=theme['secondaryColor'],
        disabled=theme['secondaryLightColor'],
        source=source,
        parent=parent,
    )
    resources.generate()

    try:
        QDir.addSearchPath('icon', resources.index)
        QDir.addSearchPath(
            'qt_material',
            os.path.join(os.path.dirname(__file__), 'resources'),
        )
    except:  # snake_case, true_property
        QDir.add_search_path('icon', resources.index)
        QDir.add_search_path(
            'qt_material',
            os.path.join(os.path.dirname(__file__), 'resources'),
        )


class ResourceGenerator:
    """"""

    # ----------------------------------------------------------------------
    def __init__(
        self,
        primary,
        secondary,
        disabled,
        source,
        parent='theme',
    ):
        """Constructor"""

        if parent.startswith('/'):
            self.index = parent
        if parent.startswith('.'):
            self.index = parent[1:]
        else:
            self.index = os.path.join(RESOURCES_PATH, parent)

        active = '#707070'

        self.contex = [
            (os.path.join(self.index, 'disabled'), disabled),
            (os.path.join(self.index, 'primary'), primary),
            (os.path.join(self.index, 'active'), active),
        ]

        self.source = source
        self.secondary = secondary

        for folder, _ in self.contex:
            shutil.rmtree(folder, ignore_errors=True)
            os.makedirs(folder, exist_ok=True)

    # ----------------------------------------------------------------------

    def generate(self):
        """"""
        for icon in os.listdir(self.source):
            if not icon.endswith('.svg'):
                continue

            with open(os.path.join(self.source, icon), 'r') as file_input:
                content_original = file_input.read()

                for folder, color in self.contex:
                    new_content = self.replace_color(content_original, color)
                    new_content = self.replace_color(
                        new_content, self.secondary, '#ff0000'
                    )

                    file_to_write = os.path.join(folder, icon)
                    with open(file_to_write, 'w') as file_output:
                        file_output.write(new_content)

    # ----------------------------------------------------------------------
    def replace_color(self, content, replace, color='#0000ff'):
        """"""
        colors = [color] + [
            ''.join(list(color)[:i] + ['\\\n'] + list(color)[i:])
            for i in range(1, 7)
        ]
        for c in colors:
            content = content.replace(c, replace)

        replace = '#ffffff00'
        color = '#000000'
        colors = [color] + [
            ''.join(list(color)[:i] + ['\\\n'] + list(color)[i:])
            for i in range(1, 7)
        ]
        for c in colors:
            content = content.replace(c, replace)

        return content
