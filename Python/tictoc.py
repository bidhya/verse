
# To measure timing
def tic():
    """
    Start the timer.
    """
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()


def toc():
    """
    Stop the timer and calculate the elapsed time.
    If tic is not called after toc, elapsed time will be cumulative from last tic call.

    Returns:
        str: A string indicating the elapsed time in seconds, minutes, or hours.
    """
    import time
    toc_str = "Toc: start time not set, call tic() first"  # default message
    if 'startTime_for_tictoc' in globals():
        elapsed_time = time.time() - startTime_for_tictoc
        if elapsed_time <= 60.0:
            toc_str = "Elapsed time is {:.3f} seconds.".format(elapsed_time)
        elif elapsed_time <= 3600.0:
            elapsed_time /= 60.0
            toc_str = "Elapsed time is {:.2f} minutes.".format(elapsed_time)
        else:
            elapsed_time /= 3600.0
            toc_str = "Elapsed time is {:.2f} hours.".format(elapsed_time)
    print(toc_str)
    return toc_str  # so we can pass it to logging etc
