#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import os
import json
from json import JSONEncoder
from pathlib import Path
from _helper_functions_ import _create_argparser_

__setting_file__ = os.path.join(str(Path.home()),
                                '.config', 'cmaq-wf', 'cmaq-wf.json')


class _Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args,
                                                                  **kwargs)
        return cls._instances[cls]


class Domain:
    def __init__(self, id, name, size, ncol, nrow):
        self.id = id
        self.type = "Domain"
        self.name = name
        self.size = size
        self.ncol = ncol
        self.nrow = nrow

    def __repr__(self):
        s = 'Domain:\n'
        s += '  Id: {}\n'.format(self.id)
        s += '  Name: {}\n'.format(self.name)
        s += '  Size: {}\n'.format(self.size)
        s += '  ncol: {}\n'.format(self.ncol)
        s += '  nrow: {}'.format(self.nrow)
        return s

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['id'], dic['name'], dic['size'],
                   dic['ncol'], dic['nrow'])


class Project:
    def __init__(self, id, name, compiler, path, path_cmaq_exe,
                 years, months, days, active=False, doms=[]):
        self.id = id
        self.type = "Project"
        self.active = active
        self.name = name
        self.compiler = compiler
        self.path = path
        self.path_cmaq_exe = path_cmaq_exe
        self.years = years
        self.months = months
        self.days = days
        self.doms = doms

    def __repr__(self):
        s = 'Name: {}\n'.format(self.name)
        s += 'Domains:\n'
        for i, d in enumerate(self.doms):
            s += '  {}- {}\n'.format(i, d.name)
        return s

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['id'], dic['name'], dic['compiler'],
                   dic['path'], dic['path_cmaq_exe'], dic['years'], dic['months'],
                   dic['days'], dic['active'], dic['doms'])

    def dom_create(self, name, size, ncol, nrow):
        ln = len(self.doms) + 1
        dom = Domain(ln, name, size, ncol, nrow)
        self.doms.append(dom)
        return dom

    def dom_delete(self, id):
        del self.doms[id]


class SettingEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


class Setting(metaclass=_Singleton):
    def __init__(self, projects=[]):
        self.type = "Setting"
        self.projects = projects

    def __repr__(self):
        s = ''
        for i, p in enumerate(self.projects):
            s += '{}- {}{}\n'.format(i, p.name,
                                     ' (active)' if p.active else '')
        return s

    @classmethod
    def defaults(cls):
        set = cls()
        set.projects = []
        proj = set.project_create(
            'cityair', 'gcc', '/mnt/disk3/projects', '/mnt/ssd2/APPS/CMAQ/',
            [2015], [1, 2, 3], list(range(1, 32)))
        proj.active = True
        proj.dom_create('eu', 36, 124, 90)
        proj.dom_create('tr', 12, 172, 94)
        proj.dom_create('aegean', 4, 103, 94)
        proj.dom_create('mediterranean', 4, 136, 97)
        proj.dom_create('central_blacksea', 4, 172, 115)
        proj.dom_create('south_central_anatolia', 4, 124, 100)
        return set

    def project_create(self, name, compiler, path, path_cmaq, years, months,
                       days):
        ln = len(self.projects) + 1
        proj = Project(ln, name, compiler, path, path_cmaq, years, months,
                       days)
        self.projects.append(proj)
        return proj

    def get_active_project(self):
        return [p for p in self.projects if p.active][0]

    @classmethod
    def load_from_file(cls, file=__setting_file__):
        try:
            with open(file, 'r') as f:
                obj = json.load(f, object_hook=cls._Decoder_)
        except FileNotFoundError:
            obj = Setting()
        projs = [p for p in obj.projects if p.active]
        if len(projs) > 1:
            for i in range(1, len(projs)):
                projs[i].active = False
        return obj

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['projects'])

    @staticmethod
    def _Decoder_(dic):
        t = dic['type']
        del dic['type']
        if t == 'Domain':
            v = Domain.fromDict(dic)
        elif t == 'Project':
            v = Project.fromDict(dic)
        elif t == 'Setting':
            v = Setting.fromDict(dic)
        return v

    def encode(self):
        return json.dumps(self, indent=2, cls=SettingEncoder)

    def decode(self):
        return json.loads(self.encode(), object_hook=self._Decoder_)

    def save(self, file=__setting_file__):
        dir = os.path.dirname(file)
        os.makedirs(dir, exist_ok=True)
        projs = [p for p in self.projects if p.active]
        if len(projs) > 1:
            for i in range(1, len(projs)):
                projs[i].active = False
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
