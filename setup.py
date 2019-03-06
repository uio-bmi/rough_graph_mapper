from setuptools import setup

setup(name='rough_graph_mapper',
      version='0.0.4',
      description='Rough Graph Mapper',
      url='http://github.com/uio-bmi/rough_graph_mapper',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      packages=["rough_graph_mapper"],
      zip_safe=False,
      install_requires=['numpy', 'mappy', 'pygssw', 'python-coveralls', 'pysam',
                        'pyfaidx', 'offsetbasedgraph', 'pyvg', 'graph_peak_caller', 'tqdm'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
        'console_scripts': ['rough_graph_mapper=rough_graph_mapper.command_line_interface:main'],
      })
