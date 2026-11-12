from importlib.metadata import PackageNotFoundError, version


def package_version() -> str:
    """Installed kolach version, or "0.0.0" when running from an uninstalled tree."""
    try:
        return version("kolach")
    except PackageNotFoundError:
        return "0.0.0"
