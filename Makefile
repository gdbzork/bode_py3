clean:
	(cd bode; rm -f *.pyc; rm -rf __pycache__)
	(cd bode/seq; rm -f *.pyc; rm -rf __pycache__)
	(cd tests; rm -f *.pyc)
	(rm -rf build)
	(rm -rf Python_3_Utility_Libraries.egg-info)
