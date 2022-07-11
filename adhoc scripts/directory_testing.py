# checking how the directory inheritance inside a function (def) works depending on the outer scope
# troubleshooting g14_linregpcr_analyze.py having problems finding file sometimes

def test_fn(x):
    return 10*x

def check_my_working_dir_inside_fn():
    
    import os
    # return 10*x
    return 'here is my working directory : ' # + os.getcwd()
    