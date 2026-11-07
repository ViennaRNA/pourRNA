if HAVE_HELP2MAN
  man1_MANS = $(ax_help2man_MANS)
  @AX_HELP2MAN_RULES@
  # generated from pourRNA --help; remove it so that make distcheck passes
  DISTCLEANFILES = $(man1_MANS)
endif
