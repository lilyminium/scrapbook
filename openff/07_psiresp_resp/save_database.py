#!/usr/bin/env python


import pathlib
import time
from typing import Optional
import tempfile
import pathlib
import glob
import shutil
import subprocess
import shutil
import os

from qcfractal.postgres_harness import find_port, FractalConfig, atexit
from qcfractal.postgres_harness import PostgresHarness as QCPostgresHarness
from qcfractal.testing import TemporaryPostgres as QCTemporaryPostgres


from qcfractal.interface.models import ResultRecord
from qcfractal.storage_sockets import storage_socket_factory
from qcfractal import FractalSnowflake, FractalSnowflakeHandler
import qcfractal.interface as ptl
import qcelemental as qcel


POSTGRES_SERVER_BACKUP = "database.psql"

class PostgresHarness(QCPostgresHarness):
    def _run(self, commands):
        command_str = " ".join(list(map(str, commands)))
        if (any(x in command_str for x in ["-p ", "--port="])
                and not any(x in command_str for x in ["-h ", "--host="])):
            commands.extend(["-h", "127.0.0.1"])

        proc = subprocess.run(commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = proc.stdout.decode()
        if not self.quiet:
            self.logger(stdout)

        ret = {"retcode": proc.returncode, "stdout": stdout, "stderr": proc.stderr.decode()}
        return ret

    def restore_database(self, filename) -> None:

        # Reasonable check here
        self._check_psql()

        self.create_database(self.config.database.database_name)

        # fmt: off
        cmds = [
            shutil.which("pg_restore"),
            "-c",
            f"--port={self.config.database.port}",
            f"--dbname={self.config.database.database_name}",
            filename
        ]
        # fmt: on

        self.logger(f"pg_backup command: {'  '.join(cmds)}")
        ret = self._run(cmds)

        if ret["retcode"] != 0:
            self.logger(ret["stderr"])
            print(ret["stderr"])
            raise ValueError("\nFailed to restore the database.\n")

class TemporaryPostgres(QCTemporaryPostgres):
    def __init__(
        self,
        database_name: Optional[str] = None,
        tmpdir: Optional[str] = None,
        quiet: bool = True,
        logger: "print" = print,
    ):
        """A PostgreSQL instance run in a temporary folder.

        ! Warning ! All data is lost when this object is deleted.

        Parameters
        ----------
        database_name : Optional[str], optional
            The database name to create.
        tmpdir : Optional[str], optional
            A directory to create the postgres instance in, if not None the data is not deleted upon shutdown.
        quiet : bool, optional
            If True, does not log any operations
        logger : print, optional
            The logger to show the operations to.
        """

        self._active = True

        if not tmpdir:
            self._db_tmpdir = tempfile.TemporaryDirectory().name
        else:
            tmpdir = pathlib.Path(tmpdir)
            self._db_tmpdir = str(tmpdir.absolute())

        self.quiet = quiet
        self.logger = logger

        config_data = {"port": find_port(), "directory": self._db_tmpdir}
        if database_name:
            config_data["database_name"] = database_name
        self.config = FractalConfig(database=config_data)
        self.psql = PostgresHarness(self.config)
        self.psql.initialize_postgres()
        self.psql.init_database()
        # self.psql.start()

        atexit.register(self.stop)


if __name__ == "__main__":
    storage = TemporaryPostgres(database_name="test_psiresp")
    storage.psql.restore_database(POSTGRES_SERVER_BACKUP)

    output_json = glob.glob("output/*.json")
    print(storage)
    records = []
    for file in output_json:
        with open(file, "r") as f:
            records.append(ResultRecord.parse_raw(f.read()))

    with FractalSnowflake(
        max_workers=1,
        storage_project_name="test_psiresp",
        storage_uri=storage.psql.database_uri(),
        reset_database=False,
        start_server=True,
    ) as tmp_server:
        print(tmp_server)
        socket = storage_socket_factory(tmp_server._storage_uri,
                                        project_name="test_psiresp",
                                        skip_version_check=True,
                                        )
        socket.add_results(records)

    storage.psql.backup_database(POSTGRES_SERVER_BACKUP)
    print(f"Saved to {POSTGRES_SERVER_BACKUP}")