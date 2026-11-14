# The report plots pick the plotly.js bundle from the
# `expressom.plotly_bundle` option (see .rep_plotly_bundle_mode()). The
# default, "auto", downloads the smaller cartesian bundle from cdn.plot.ly the
# first time a session builds a plot -- the test suite must never depend on
# network access, so pin the full (locally installed) bundle here. The
# download / swap machinery itself is unit-tested with a local dummy file in
# test-report-load-speed.R.
options(expressom.plotly_bundle = "full")
