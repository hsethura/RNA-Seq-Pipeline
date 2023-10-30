def format_table(lst_rows, column_sep = " ", right_justify=None):
    column_widths = [0] * len(lst_rows[0])
    for row in lst_rows:
        for i in xrange(len(row)): column_widths[i] = max([column_widths[i], len(row[i])])
    if right_justify is None:
        right_justify = [False] * len(column_widths)
    justify = ['>' if j else '' for j in right_justify]
    row_format = column_sep.join(['{:' + justify + str(width) + 's}'  for (width, justify) in zip(column_widths, justify)])
    lst_lines = [row_format.format(*lst_rows[0])]
    header_sep = '-' * (sum(column_widths) + len(column_sep) * (len(column_widths) - 1))
    lst_lines.append(header_sep)
    lst_lines.extend([row_format.format(*row) for row in lst_rows[1:]])
    return lst_lines

