from amuse.support.literature import TrackLiteratureReferences, LiteratureReferencesMixIn

class OMUSE(LiteratureReferencesMixIn):
    """
        .. [#] OMUSE 1.0: a framework for multi-model oceanographic simulations, Pelupessy et al., Geoscientific Model Development, 10, 3167
    """

TrackLiteratureReferences.default().register_class(OMUSE)
