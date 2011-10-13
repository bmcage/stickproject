#
# Copyright (C) 2005-2007  Donald N. Allingham
# Copyright (C) 2008-2009  Gary Burton 
# Copyright (C) 2009       Doug Blank <doug.blank@gmail.com>
# Copyright (C) 2009       Benny Malengier <bm@cage.ugent.be>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


"""
This package implements access to configuration.
"""

#---------------------------------------------------------------
#
# System imports
#
#---------------------------------------------------------------
import os
import time
import ConfigParser
import errno

#---------------------------------------------------------------
#
# local imports
#
#---------------------------------------------------------------


#---------------------------------------------------------------
#
# Constants
#
#---------------------------------------------------------------


#---------------------------------------------------------------
#
# Local functions
#
#---------------------------------------------------------------
def eval_item(setting):
    """
    Given a value from an ini file, return it in its proper type.
    May be recursively called, in the case of nested structures.
    """
    setting = setting.strip()
    value = None
    if setting.startswith("'") and setting.endswith("'"):
        value = setting[1:-1]
    elif setting.startswith("[") and setting.endswith("]"):
        list_data = setting[1:-1].replace('\n', '').replace(' ','')
        if setting.startswith('[['):
            # we do not allow string data, list of list of numbers
            value = eval(setting)
        else:
            if list_data:
                value = [eval_item(item) for item in list_data.split(",")]
            else:
                value = []
    elif setting == "True":
        value = True 
    elif setting == "False":
        value = False
    elif "." in setting or "e" in setting or "E" in setting:
        value = float(setting)
    else:
        value = int(setting)
    return value

#---------------------------------------------------------------
#
# Classes
#
#---------------------------------------------------------------
class ConfigManager(object):
    """
    Abstract base class to construct a singleton CONFIGMAN where all 
    settings are stored.
    """
    
    def __init__(self, filename):
        """ 
        Configure manager constructor takes an optional filename.

        The data dictionary stores the settings:

           self.data[section][setting] = value

        The value has a type that matches the default. It is an error
        to attempt to set the setting to a different type. To change
        the type, you must re-register the setting, and re-set the
        value.

        The default values are given in Python code and stored here
        on start-up:

           self.default[section][setting] = default_value

        Callbacks are stored as callables here:

           self.callbacks[section][setting] = (id, func)

        The default filename (usually the one you are reading from)
        is stored as self.filename. However, you can save to another
        filename using self.save(otherfilename).
        """
        self._cb_id = 0 # callback id counter
        self.filename = filename
        self.callbacks = {}
        self.default = {}
        self.help = {}
        self.data = {}
        self.reset()
        self.register_defaults()
        self.load()

    def register_defaults(self):
        """
        Must be implemented in inheriting classes
        """
        raise NotImplementedError

    def reset(self, key=None):
        """
        Resets one, a section, or all settings values to their defaults.
        This does not disconnect callbacks.
        """
        if key is None:
            section = None
            setting = None
        elif "." in key:
            section, setting = key.split(".", 1)
        else: # key is not None and doesn't have a "."
            section = key
            setting = None
        # Now, do the reset on the right parts:
        if section is None:
            self.data = {}
            for section in self.default:
                self.data[section] = {}
                for setting in self.default[section]:
                    self.data[section][setting] = self.default[section][setting]
        elif setting is None:
            self.data[section] = {}
            for setting in self.default[section]:
                self.data[section][setting] = self.default[section][setting]
        else:
            self.data[section][setting] = self.default[section][setting]
        # Callbacks are still connected

    def get_sections(self):
        """
        Return all section names.
        """
        return self.data.keys()

    def get_section_settings(self, section):
        """
        Return all section setting names.
        """
        return self.data[section].keys()

    def load(self):
        """ 
        Loads an .ini into self.data.
        """
        if self.filename and os.path.exists(self.filename):
            parser = ConfigParser.ConfigParser()
            parser.read(self.filename)
            for sec in parser.sections():
                name = sec.lower()
                if name not in self.data:
                    # Add the setting from file
                    self.data[name] = {}
                for opt in parser.options(sec):
                    setting = parser.get(sec, opt).strip()
                    value = eval_item(setting)
                    #Now, let's test and set:
                    if name in self.default and opt.lower() in self.default[name]:
                        if type(value) == type(self.default[name][opt.lower()]):
                            self.data[name][opt.lower()] = value
                        else:
                            print ("WARNING: ignoring key with wrong type "
                                   "'%s.%s'" % (name, opt.lower()))
                    else:
                        self.data[name][opt.lower()] = value
        else:
            raise ValueError, 'The file %s cannot be found' % self.filename

    def save(self, filename=None, helpafter=False):
        """
        Saves the current section/settings to an .ini file. Optional filename
        will override the default filename to save to, if given.
        help is by default under the option, if helpafter=True, it is behind.
        """
        if filename is None:
            filename = const.INIFILE_DEFAULT
        if filename:
            try:
                head = os.path.split( filename )[0]
                os.makedirs( head )
            except OSError, exp:
                if exp.errno != errno.EEXIST:
                    raise
            key_file = open(filename, "w")
            key_file.write(";; StickProject key file\n")
            key_file.write((";; Automatically created at %s" % 
                      time.strftime("%Y/%m/%d %H:%M:%S")) + "\n\n")
            sections = sorted(self.data)
            for section in sections:
                key_file.write(("[%s]\n") % section)
                keys = sorted(self.data[section])
                for key in keys:
                    help = ""
                    value = self.data[section][key]
                    default = ""
                    if self.has_default("%s.%s" % (section, key)):
                        help = self.help[section][key]
                        if value == self.get_default("%s.%s" % (section, key)):
                            default = ";;"
                        else:
                            default = ""
                    if isinstance(value, long):
                        value = int(value)
                    if (help and helpafter) or default:
                        key_file.write(("%s%s=%s")% (default, key, repr(value)))
                    else:
                        key_file.write(("%s%s=%s\n")% (default, key, repr(value)))
                    
                    if help and default:
                        key_file.write(" ;; %s\n" % help)
                    elif help:
                        key_file.write(";; %s\n" % help)
                    elif default:
                        key_file.write('\n')
                key_file.write("\n")
            key_file.close()
        # else, no filename given; nothing to save so do nothing quietly

    def get(self, key):
        """
        Get the setting's value. raise an error if an invalid section.setting.
        Key is a sting in the "section.setting" format.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            raise AttributeError("Invalid config section.setting name: '%s'" % 
                                 key)
        if section not in self.data:
            raise AttributeError("No such config section name: '%s'" % section)
        if setting not in self.data[section]:
            raise AttributeError("No such config setting name: '%s.%s'" % 
                                 (section, setting))
        return self.data[section][setting]

    def is_set(self, key):
        """
        Does the setting exist? Returns True if does, False otherwise.
        Key is a sting in the "section.setting" format.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            return False
        if section not in self.data:
            return False
        if setting not in self.data[section]:
            return False
        return True

    def get_or_fallback(self, key, fallback):
        """
        get a key, if not present, return fallback instead
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            return fallback
        if section not in self.data:
            return fallback
        if setting not in self.data[section]:
            return fallback
        return self.data[section][setting]

    def has_default(self, key):
        """
        Does the setting have a default value? Returns True if it does, 
        False otherwise. Key is a sting in the "section.setting" format.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            return False
        if section not in self.default:
            return False
        if setting not in self.default[section]:
            return False
        return True

    def get_default(self, key):
        """
        Get the setting's default value. Raises an error if invalid key is
        give. Key is a sting in the "section.setting" format.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            raise AttributeError("Invalid config section.setting name: '%s'" % 
                                 key)
        if section not in self.default:
            raise AttributeError("No such config section name: '%s'" % section)
        if setting not in self.default[section]:
            raise AttributeError("No such config setting name: '%s.%s'" % 
                                 (section, setting))
        return self.default[section][setting]

    def register(self, key, default, help=''):
        """
        Register a section.setting, and set the default.
        help is an optional help text
        Will overwrite any previously set default, and set setting if not one.
        The default value deterimines the type of the setting.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            raise AttributeError("Invalid config section.setting name: '%s'" % 
                                 key)
        if section not in self.data:
            self.data[section] = {}
        if section not in self.default:
            self.default[section] = {}
        if section not in self.help:
            self.help[section] = {}
        if section not in self.callbacks:
            self.callbacks[section] = {}
        if setting not in self.callbacks[section]:
            self.callbacks[section][setting] = []
        # Add the default value to settings, if not exist:
        if setting not in self.data[section]:
            self.data[section][setting] = default
        # Set the default, regardless:
        self.default[section][setting] = default
        # Set the help text
        self.help[section][setting] = help

    def connect(self, key, func):
        """
        Connect a callback func that gets called when key is changed.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            raise AttributeError("Invalid config section.setting name: '%s'" % 
                                 key)
        if section not in self.data:
            raise AttributeError("No such config section name: '%s'" % section)
        if setting not in self.data[section]:
            raise AttributeError("No such config setting name: '%s.%s'" % 
                                 (section, setting))
        self._cb_id += 1
        self.callbacks[section][setting].append((self._cb_id, func))
        return self._cb_id

    def disconnect(self, callback_id):
        """
        Removes a callback given its callback ID. The ID is generated and
        returned when the function is connected to the key (section.setting).
        """
        for section in self.callbacks:
            for setting in self.callbacks[section]:
                for (cbid, func) in self.callbacks[section][setting]:
                    if callback_id == cbid:
                        self.callbacks[section][setting].remove((cbid, func))

    def emit(self, key):
        """
        Emits the signal "key" which will call the callbacks associated
        with that setting.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            raise AttributeError("Invalid config section.setting name: '%s'" % 
                                 key)
        if section not in self.callbacks:
            raise AttributeError("No such config section name: '%s'" % section)
        if setting not in self.callbacks[section]:
            raise AttributeError("No such config setting name: '%s.%s'" % 
                                 (section, setting))
        for (cbid, func) in self.callbacks[section][setting]:
            func(self, 0, str(self.data[section][setting]), None) 

    def set(self, key, value):
        """
        Set the setting's value. There are only two ways to get into
        the data dictionary: via the load() method that reads a file,
        or from this method.
        """
        if "." in key:
            section, setting = key.split(".", 1)
        else:
            raise AttributeError("Invalid config section.setting name: '%s'" % 
                                 key)
        if section not in self.data:
            raise AttributeError("No such config section name: '%s'" % section)
        if setting not in self.data[section]:
            raise AttributeError("No such config setting name: '%s.%s'" % 
                                 (section, setting))
        # Check value to see if right type:
        if type(value) == long:
            value = int(value)
        if type(value) == unicode:
            value = str(value)
        if self.has_default(key):
            if type(self.get_default(key)) != type(value):
                raise AttributeError("attempting to set '%s' to wrong type "
                                     "'%s'; should be '%s'" %
                                     (key, type(value), 
                                      type(self.get_default(key))))
        if (setting in self.data[section] and 
            self.data[section][setting] == value):
            # Do nothing if existed and is the same
            pass
        else:
            # Set the value:
            self.data[section][setting] = value
            # Only call callback if the value changed!
            if (section in self.callbacks and 
                setting in self.callbacks[section]):
                for (cbid, func) in self.callbacks[section][setting]:
                    # Call all callbacks for this key:
                    func(self, 0, str(self.data[section][setting]), None) 
