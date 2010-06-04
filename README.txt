Installation
============

Before installing, make a copy of *org_settings.cfg.sample* to *org_settings.cfg* :

  $> cp org_settings.cfg.sample org_settings.cfg

In the new *org_settings.cfg*, create/edit the paths and categories desired for
your system as appropriate.  When you have configured the file to your
satisfaction, copy it into the root source directory:

  $> cp org_settings.cfg src/chipsequtil/

You can then install the package with:

  $> python setup.py install


If you'd like to install the package to a non-system directory (e.g., if you
don't have permission to install system-wide packages), you can provide the
*--prefix=PATH* argument to the install command:

  $> python setup.py install --prefix=/path/to/dir

Remember to add */path/to/dir* to your PYTHONPATH environment variable if it
is not already there.  If you wish to add more system-wide paths/organisms to
org_settings.cfg, either edit the file in the source directory as above and
reinstall (good way) or edit the file in the directory where the package is
installed (less good way).
