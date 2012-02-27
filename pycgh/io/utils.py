def _file_handle(file_ref, mode='r'):

    if not mode in 'rw':
        raise ValueError("mode must be 'r' or 'w'")

    def _is_string_like(obj):
        try:
            obj + ''
        except (TypeError, ValueError):
            return False
        return True

    try:
        if _is_string_like(file_ref):
            if file_ref.endswith('.gz'):
                import gzip
                fh = gzip.open(file_ref, mode='%sb' % mode)
            else:
                if mode == 'r':
                    fh = open(file_ref, 'U')
                else:
                    fh = open(file_ref, 'w')
        else:
            fh = file_ref
    except TypeError:
        raise ValueError('input file must be a path or file handle')

    return fh
