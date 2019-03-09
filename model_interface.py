""""
Definition of model interface
"""


class IModel(object):

    def name(self) -> str:
        pass

    def create_model(self, active_descriptors: str, inactive_descriptors: str, model_name: str) -> dict:
        raise NotImplemented()

    def save_to_json_file(self, output_file: str, model: dict):
        raise NotImplemented()

    def score_model(self, data: list, fragments_file: str, descriptors_file: str, output_file: str):
        raise NotImplemented()


