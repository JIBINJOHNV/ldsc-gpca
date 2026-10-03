"""Standard-library-only help formatting, shared with R's argparse wrapper."""
import argparse
import re


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    def __init__(self, prog):
        super().__init__(prog, width=96, max_help_position=32)

    def _get_help_string(self, action):
        text = action.help or ''
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)) or '--help' in action.option_strings:
            return text
        if action.required:
            if 'required; no default' not in text.lower():
                text += ' [Required; no default.]'
        elif action.option_strings and not re.search(r'\bdefault\s*:', text, re.I):
            default = ('unset' if action.default is None else
                       'on' if action.default is True else
                       'off' if action.default is False else '%(default)s')
            text += f' [Default: {default}.]'
        if action.choices:
            text += ' [Choices: ' + ', '.join(map(str, action.choices)) + '.]'
        return text

    def _metavar_formatter(self, action, default_metavar):
        if action.choices and action.option_strings and action.metavar is None:
            return lambda size: ('CHOICE',) * size
        return super()._metavar_formatter(action, default_metavar)

    def _get_default_metavar_for_optional(self, action):
        return {int: 'INTEGER', float: 'FLOAT'}.get(action.type, action.dest.upper())

    def _fill_text(self, text, width, indent):
        # Wrap prose, while preserving multiline schemas and shell examples.
        if '\n' not in text:
            return argparse.HelpFormatter._fill_text(self, text, width, indent)
        return super()._fill_text(text, width, indent)

    def _format_action(self, action):
        return super()._format_action(action) + '\n'
