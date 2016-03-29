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
        import sys
        if not os.path.isfile(os.path.join(os.path.split(__file__)[0], 'blackbox', 'data')):
            print 'recording git commit'
            from blackbox.accessoryFunctions import make_path
            make_path(os.path.join(os.path.split(__file__)[0], 'blackbox', 'data'))
            with open(os.path.join(os.path.split(__file__)[0], 'blackbox', 'data', 'git.dat'), 'w') as git:
                git.write(os.popen('git rev-parse --short HEAD').read().rstrip())
        # Attempt to detect whether we were called from setup() or by another
        # command.  If we were called by setup(), our caller will be the
        # 'run_command' method in 'distutils.dist', and *its* caller will be
        # the 'run_commands' method.  If we were called any other way, our
        # immediate caller *might* be 'run_command', but it won't have been
        # called by 'run_commands'.  This is slightly kludgy, but seems to
        # work.
        #
        caller = sys._getframe(2)
        caller_module = caller.f_globals.get('__name__', '')
        caller_name = caller.f_code.co_name
        if caller_module != 'distutils.dist' or caller_name != 'run_commands':
            # We weren't called from the command line or setup(), so we
            # should run in backward-compatibility mode to support bdist_*
            # commands.
            install.run(self)
        else:
            self.do_egg_install()


setup(
    name='docker-assembly',
    version='0.0.dev1',
    packages=['blackbox'],
    url='https://github.com/MikeKnowles/docker-assembly',
    package_data=dict(MBBspades=['blackbox/data/*.dat']),
    include_package_data=True,
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='Assembly pipeline using Docker and SPAdes',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0',
                      'pyyaml'],
    scripts=['bin/MBBspades'],
    cmdclass=dict(install=RecordGit)
)
