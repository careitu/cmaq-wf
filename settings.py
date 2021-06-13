#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import os
import json
import warnings
from json import JSONEncoder
from pathlib import Path
from _helper_functions_ import _create_argparser_

__setting_file__ = os.path.join(str(Path.home()),
                                '.config', 'cwf', 'cwf.json')


class _Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args,
                                                                  **kwargs)
        return cls._instances[cls]


class Domain:
    def __init__(self, id, id2, name, size, ncol, nrow, sub_doms=None):
        self.id = id
        self.id2 = id2
        self.type = 'Domain'
        self.name = name
        self.size = size
        self.ncol = ncol
        self.nrow = nrow
        self.__parent__ = None
        self.doms = {} if sub_doms is None else sub_doms

    def __repr__(self):
        s = 'Domain:\n'
        s += '  Id: {} ({})\n'.format(self.id, self.id2)
        s += '  Name: {}\n'.format(self.name)
        s += '  Size: {}\n'.format(self.size)
        s += '  ncol: {}\n'.format(self.ncol)
        s += '  nrow: {}\n'.format(self.nrow)
        if len(self.doms) > 0:
            s += '  Sub-Domains:\n'
            for k, v in self.doms.items():
                s += '    {}- {}\n'.format(v.id, k)
        return s

    def get_dom_by_id(self, id):
        dom = None
        if self.id == id:
            dom = self
        if dom is None and len(self.doms) > 0:
            for v in self.doms.values():
                dom = v.get_dom_by_id(id)
                if dom is not None:
                    break
        return dom

    def get_sub_ids(self):
        ids = [self.id]
        if len(self.doms) > 0:
            for v in self.doms.values():
                ids = ids + v.get_sub_ids()
        return ids

    def append(self, id2, name, size, ncol, nrow):
        ln = len(self.doms) + 1
        dom = Domain(int(f"{self.id}{ln}"), id2, name, size, ncol, nrow)
        dom.__parent__ = self
        self.doms[name] = dom
        return dom

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['id'], dic['id2'], dic['name'], dic['size'],
                   dic['ncol'], dic['nrow'], dic['doms'])


class Project:
    def __init__(self, id, name, compiler, path, path_cmaq_exe,
                 years, months, days, active=False, doms=None):
        self.id = id
        self.type = 'Project'
        self.active = active
        self.name = name
        self.compiler = compiler
        self.path = path
        self.path_cmaq_exe = path_cmaq_exe
        self.years = years
        self.months = months
        self.days = days
        self.doms = doms
        self.doms = {} if doms is None else doms

    def __repr__(self):
        s = 'Name: {}\n'.format(self.name)
        s += 'Domains:\n'
        for k, v in self.doms.items():
            s += '  {}- {}\n'.format(v.id, k)
        return s

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['id'], dic['name'], dic['compiler'],
                   dic['path'], dic['path_cmaq_exe'], dic['years'],
                   dic['months'], dic['days'], dic['active'], dic['doms'])

    def get_dom_by_id(self, id):
        dom = None
        if dom is None and len(self.doms) > 0:
            for v in self.doms.values():
                dom = v.get_dom_by_id(id)
                if dom is not None:
                    break
        return dom

    def get_dom_ids(self):
        ids = []
        for v in self.doms.values():
            ids = ids + v.get_sub_ids()
        return ids

    def append(self, id2, name, size, ncol, nrow):
        ln = len(self.doms) + 1
        dom = Domain(ln, id2, name, size, ncol, nrow)
        self.doms[name] = dom
        return dom

    def delete_dom(self, id):
        del self.doms[id]


class SettingEncoder(JSONEncoder):
    def default(self, o):
        dic = {k: v for k, v in o.__dict__.items() if not k.startswith('__')}
        return dic


class Setting(metaclass=_Singleton):
    def __init__(self, projects={}):
        self.type = 'Setting'
        self.projects = projects

    def __repr__(self):
        s = ''
        for k, v in self.projects.items():
            s += '{}- {}{}\n'.format(v.id, k,
                                     ' (active)' if v.active else '')
        return s

    @classmethod
    def defaults(cls):
        set = cls()
        proj = set.create_proj(
            'cityair', 'gcc', '/mnt/disk3/projects', '/mnt/ssd2/APPS/CMAQ/',
            [2015], [1, 2, 3], list(range(1, 32)))
        proj.active = True
        eu = proj.append(1, 'eu', 36, 124, 90)
        tr = eu.append(2, 'tr', 12, 172, 90)
        tr.append(3, 'aegean', 4, 103, 94)
        tr.append(4, 'mediterranean', 4, 136, 97)
        tr.append(5, 'central_blacksea', 4, 172, 115)
        tr.append(6, 'south_central_anatolia', 4, 124, 100)
        return set

    def create_proj(self, name, compiler, path, path_cmaq, years, months,
                    days):
        ln = len(self.projects) + 1
        proj = Project(ln, name, compiler, path, path_cmaq, years, months,
                       days)
        self.projects[name] = proj
        return proj

    def get_active_proj(self, warn=True):
        keys = list(self.projects.keys())
        projs = {k: v for k, v in self.projects.items() if v.active}
        if warn:
            if len(projs) > 1:
                msg = "You have multiple active projects. Returning first."
                warnings.warn(msg)
        if len(projs) == 0:
            warnings.warn("You don't have any project")
        return projs[keys[0]]

    @staticmethod
    def set_parents(dom):
        if len(dom.doms) > 0:
            for k, v in dom.doms.items():
                v.__parent__ = dom
                Setting.set_parents(v)

    @classmethod
    def load_from_file(cls, file=__setting_file__):
        try:
            with open(file, 'r') as f:
                obj = json.load(f, object_hook=cls._Decoder_)
        except FileNotFoundError:
            obj = Setting()
        projs = {k: v for k, v in obj.projects.items() if v.active}
        if len(projs) > 1:
            keys = list(projs.keys())
            for i in range(1, len(keys)):
                projs[keys[i]].active = False
        # set parents
        for k, v in obj.projects.items():
            for k2, v2 in v.doms.items():
                Setting.set_parents(v2)

        return obj

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['projects'])

    @staticmethod
    def _Decoder_(dic):
        if 'type' in dic.keys():
            t = dic['type']
            del dic['type']
            if t == 'Domain':
                v = Domain.fromDict(dic)
            elif t == 'Project':
                v = Project.fromDict(dic)
            elif t == 'Setting':
                v = Setting.fromDict(dic)
            return v
        return dic

    def encode(self):
        return json.dumps(self, indent=2, cls=SettingEncoder)

    def decode(self):
        return json.loads(self.encode(), object_hook=self._Decoder_)

    def save(self, file=__setting_file__):
        dir = os.path.dirname(file)
        os.makedirs(dir, exist_ok=True)
        projs = {k: v for k, v in self.projects.items() if v.active}
        if len(projs) > 1:
            keys = list(projs.keys())
            for i in range(1, len(keys)):
                projs[keys[i]].active = False
        with open(file, 'w') as f:
            json.dump(self, f, indent=2, cls=SettingEncoder)

    def load(self, file=__setting_file__):
        self = self.load_from_file(file)


setting = Setting.load_from_file()

if __name__ == "__main__":
    DESCRIPTION = 'Save/load CMAQ Settings\n'
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -d \n'

    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-d', '--default', help="create default settings file",
                   default=False, action="store_true")
    args = p.parse_args()
    if args.default:
        os.remove(__setting_file__)
        Setting.defaults().save()
        print('Default settings were saved at {}'.format(__setting_file__))
    else:
        print(setting.encode())
