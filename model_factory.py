import model_interface

_required_models = {}


def register_model(name, constructor):
    _required_models[name] = constructor


def create_model(name: str) -> model_interface.IModel:
    if _required_models == {}:
        _discover_models()
    return _required_models[name]


def _discover_models():
    import index_model
    import descriptors_model

