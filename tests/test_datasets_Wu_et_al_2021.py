import logging

logging.basicConfig(level=logging.DEBUG, format="%(levelname)s | %(asctime)s | %(message)s")
import unittest
from pyprism.datasets import WuEtAl2021, wu_et_al_2021_processed
from pyprism.store import Store, GSM
from pathlib import Path
from tempfile import TemporaryDirectory
import os


class WuEtAl2021(unittest.TestCase):

    def test_full_download(self):
        with TemporaryDirectory() as tmpdirname:
            store = Store(path=Path(tmpdirname))
            WuEtAl2021(store=store).get()
            # for filename in os.listdir(tmpdirname):
            #     filepath = store.get_path().joinpath(filename)
            #     self.assertEqual(sha256sum(filepath), filename[:64])

    def test_partial_download(self):
        with TemporaryDirectory() as tmpdirname:
            store = Store(path=Path(tmpdirname))
            logging.debug(f"Store: {store}")
            data = WuEtAl2021(store=store).get(gsm=GSM("GSM5354522"))
            logging.info(data)
            logging.debug(f" All files in store: {os.listdir(tmpdirname)}")


class WuEtAl2021Processed(unittest.TestCase):

    def test_full_run(self):
        with TemporaryDirectory() as tmpdirname:
            store = Store(path=Path(tmpdirname))
            store_element = wu_et_al_2021_processed(store=store)
            dataset = store_element.get()
            logging.info(dataset)


if __name__ == '__main__':
    unittest.main()
