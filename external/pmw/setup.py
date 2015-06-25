#!/usr/bin/env python
# setup.py 
from distutils.core import setup 


setup(name="Pmw",
      version='1.3.2',
      description = 'Python Mega Widgets',
      author="Telstra Corporation Limited, Australia",
      author_email="",
      url='http://pmw.sourceforge.net/',
      
      package_dir = { "Pmw":"Pmw"},
      
      packages=['Pmw', 'Pmw.Pmw_1_3',
		'Pmw.Pmw_1_3.lib',],

      package_data={'Pmw': ['Pmw_1_3/lib/Pmw.def',
                            'Pmw_1_3/doc/*.html',
                            'Pmw_1_3/doc/*.gif',
                            'Pmw_1_3/doc/*.py',
                            'Pmw_1_3/contrib/*.py',
                            'Pmw_1_3/contrib/README',
                            'Pmw_1_3/demos/*.py',
                            'Pmw_1_3/tests/*.py',
                            'Pmw_1_3/tests/*.gif',
                            'Pmw_1_3/tests/*.bmp',
                            'Pmw_1_3/bin/*.py',
			   ]
                   },
      
      license='BSD',
      long_description='''Pmw is a toolkit for building high-level compound widgets, or megawidgets, 
	constructed using other widgets as component parts. It promotes consistent look and feel within
	 and between graphical applications, is highly configurable to your needs and is easy to use.''',
      classifiers = [
          'Development Status :: Alpha',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'Intended Audience :: System Administrators',
          'License :: BSD',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: GUI',
          ],
     )
