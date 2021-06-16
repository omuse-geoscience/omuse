from amuse.support.literature import TrackLiteratureReferences, LiteratureReferencesMixIn

class OMUSE(LiteratureReferencesMixIn):
    """
        .. [#] OMUSE 1.0: a framework for multi-model oceanographic simulations, Pelupessy et al., Geoscientific Model Development, 10, 3167
    """

    @classmethod
    def version(cls):
        """
        Overide the class-method from AMUSE, so we don't accidentally report
        the AMUSE version instead of the OMUSE version.
        """

        try:
            from omuse.version import version
        except (ImportError):
            version = "unknown"

        return version

TrackLiteratureReferences.default().register_class(OMUSE)
