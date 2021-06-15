#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import os as _os
import json as _json
from warnings import warn as _warn
from json import JSONEncoder as _je
from pathlib import Path as _path

__setting_file__ = _os.path.join(str(_path.home()),
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


class Paths():
    def __init__(self, **kwargs):
        self.type = 'Paths'
        self.__dict__.update(kwargs)

    def __repr__(self):
        s = 'Paths:\n'
        for k, v in self.__dict__.items():
            s += '  {}: {}\n'.format(k, v)
        return s

    @classmethod
    def fromDict(cls, dic):
        return cls(**dic)


class Project:
    def __init__(self, id, name, compiler, cmaq_ver, years, months, days,
                 paths=None, active=False, doms=None):
        self.id = id
        self.type = 'Project'
        self.active = active
        self.name = name
        self.compiler = compiler
        self.cmaq_ver = cmaq_ver
        self.path = paths
        self.years = years
        self.months = months
        self.days = days
        self.doms = doms
        self.doms = {} if doms is None else doms
        self.__setting__ = None

    def __repr__(self):
        s = 'Name: {}\n'.format(self.name)
        s += 'Domains:\n'
        for k, v in self.doms.items():
            s += '  {}- {}\n'.format(v.id, k)
        return s

    @classmethod
    def fromDict(cls, dic):
        return cls(dic['id'], dic['name'], dic['compiler'], dic['cmaq_ver'],
                   dic['years'], dic['months'], dic['days'], dic['path'],
                   dic['active'], dic['doms'])

    def activate(self):
        self.__setting__.activate(self.name)

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


class SettingEncoder(_je):
    def default(self, o):
        dic = {k: v for k, v in o.__dict__.items() if not k.startswith('__')}
        return dic


class Setting(metaclass=_Singleton):
    def __init__(self, projects=None):
        self.type = 'Setting'
        self.projects = {} if projects is None else projects

    def __repr__(self):
        s = ''
        for k, v in self.projects.items():
            s += '{}- {}{}\n'.format(v.id, k,
                                     ' (active)' if v.active else '')
        return s

    @classmethod
    def defaults(cls):
        from os.path import join as _join
        from copy import deepcopy as _dcp
        set = cls()
        set.projects = {}
        proj_name = 'cityair'
        dir_proj = _join('/mnt/disk3/projects/', proj_name)
        dir_cmaq_app = '/mnt/ssd2/APPS/CMAQ'
        proj = set.new_proj('cityair', 'gcc', '532', [2015], [1, 2, 3],
                            list(range(1, 32)), active=False)
        proj.path = Paths(proj=dir_proj,
                          cmaq_app=dir_cmaq_app,
                          wps=_join(dir_proj, 'WPS'),
                          wrf=_join(dir_proj, 'WRF'),
                          cctm=_join(dir_proj, 'cmaq'),
                          mcip=_join(dir_proj, 'mcip'),
                          icon=_join(dir_proj, 'icon'),
                          bcon=_join(dir_proj, 'bcon'),
                          emis=_join(dir_proj, 'emis'))
        eu = proj.append(1, 'eu', 36, 124, 90)
        tr = eu.append(2, 'tr', 12, 172, 90)
        tr.append(3, 'aegean', 4, 103, 94)
        tr.append(4, 'mediterranean', 4, 136, 97)
        tr.append(5, 'central_blacksea', 4, 172, 115)
        tr.append(6, 'south_central_anatolia', 4, 124, 100)
        #
        proj_name = 'test_cityair'
        dir_proj = _join('/mnt/disk3/projects/', proj_name)
        proj2 = _dcp(proj)
        proj2.active = True
        proj2.id = 2
        proj2.name = 'test_cityair'
        proj2.path = Paths(proj=dir_proj,
                           cmaq_app=dir_cmaq_app,
                           wps=_join(dir_proj, 'WPS'),
                           wrf=_join(dir_proj, 'WRF'),
                           cctm=_join(dir_proj, 'cmaq'),
                           mcip=_join(dir_proj, 'mcip'),
                           icon=_join(dir_proj, 'icon'),
                           bcon=_join(dir_proj, 'bcon'),
                           emis=_join(dir_proj, 'emis'))
        set.projects[proj_name] = proj2
        return set

    def append(self, proj):
        proj.id = len(self.projects) + 1
        self.projects[proj.name] = proj

    def new_proj(self, name, compiler, cmaq_ver, years,
                 months, days, paths=None, active=False, doms=None):
        proj = Project(None, name, compiler, cmaq_ver, years, months, days,
                       paths, active, doms)
        proj.__setting__ = self
        self.append(proj)
        return proj

    def activate(self, proj_name):
        keys = list(self.projects.keys())
        if proj_name in keys:
            for k, v in self.projects.items():
                v.active = k == proj_name
            self.save()
        else:
            raise ValueError('{} is not a project name'.format(proj_name))

    def get_active_proj(self, warn=True):
        projs = {k: v for k, v in self.projects.items() if v.active}
        keys = list(projs.keys())
        if warn:
            if len(projs) > 1:
                msg = "You have multiple active projects. Returning first."
                _warn(msg)
        if len(projs) == 0:
            raise ValueError("You don't have any active project")
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
                obj = _json.load(f, object_hook=cls._Decoder_)
        except FileNotFoundError:
            obj = Setting()
        except _json.decoder.JSONDecodeError:
            print('Settings cannot be loaded.')
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
            elif t == 'Paths':
                v = Paths.fromDict(dic)
            return v
        return dic

    def encode(self):
        return _json.dumps(self, indent=2, cls=SettingEncoder)

    def decode(self):
        return _json.loads(self.encode(), object_hook=self._Decoder_)

    def save(self, file=__setting_file__):
        dir = _os.path.dirname(file)
        _os.makedirs(dir, exist_ok=True)
        projs = {k: v for k, v in self.projects.items() if v.active}
        if len(projs) > 1:
            keys = list(projs.keys())
            for i in range(1, len(keys)):
                projs[keys[i]].active = False
        with open(file, 'w') as f:
            _json.dump(self, f, indent=2, cls=SettingEncoder)

    def load(self, file=__setting_file__):
        self = self.load_from_file(file)


setting = Setting.load_from_file()

if __name__ == "__main__":
    from _helper_functions_ import _create_argparser_

    DESCRIPTION = 'Save/load CMAQ Settings\n'
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -d \n'

    p = _create_argparser_(DESCRIPTION, EPILOG)
    g = p.add_mutually_exclusive_group(required=False)
    g.add_argument('-d', '--default', help="create default settings",
                   default=False, action="store_true")
    g.add_argument('-a', '--activate', metavar='PROJECT_NAME',
                   help="Activate a project")
    g.add_argument('-c', '--create', help="create folders for active project",
                   default=False, action="store_true")

    args = p.parse_args()

    if args.activate is not None:
        setting.activate(args.activate)
        print('{} was activated'.format(args.activate))

    if args.create:
        proj = setting.get_active_proj()
        for p in proj.path.__dict__.values():
            _os.makedirs(p, exist_ok=True)
            print("* '{}' created.".format(p))

    if args.default:
        if _os.path.isfile(__setting_file__):
            _os.remove(__setting_file__)
        Setting.defaults().save()
        print('Default settings were saved at {}'.format(__setting_file__))

    if args.print:
        print(setting.encode())

    if not args.create and not args.print and \
       args.activate is None and not args.default:
        print('CWF Projects:\n{}'.format(setting))