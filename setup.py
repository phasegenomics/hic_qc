import versioneer
from distutils.core import setup

setup(
    name='hic_qc',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='GNU Affero General Public License v3.0'
)
