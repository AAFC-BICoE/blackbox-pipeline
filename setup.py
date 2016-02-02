try:
    from setuptools import setup, Command
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup, Command
    from distutils.command.install import install
import os


class RecordGit(install):
    description = 'include git SHA-1 sum'

    def run(self):
        print 'recording git commit'
        with open(os.path.join(os.path.split(__file__)[0], 'MBBSpades', 'data', 'git.dat')) as git:
            git.write(os.popen('git rev-parse --short HEAD').read().rstrip())
        install.run(self)

setup(
    name='docker-assembly',
    version='0.0dev1',
    packages=['MBBspades'],
    url='https://github.com/MikeKnowles/docker-assembly',
    package_data=dict(MBBspades=['MBBspades/data/*.dat']),
    include_package_data=True,
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='Assembly pipeline using Docker and SPAdes',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0'],
    cmdclass=dict(install=RecordGit)
)
