from setuptools import setup

requires = [
    'pylab',
    'numpy'
    ]

setup(name='phhotoreceptor',
      version='0.1',
      description='A HH insect photoreceptor simulator',
      url='',
      author='Francisco J. Hernandez Heras',
      author_email='fjhheras@gmail.com',
      license='GPLv3',
      packages=['phhotoreceptor'],
      install_requires=[
          'matplotlib',
          'numpy'
      ],
      zip_safe=False)
