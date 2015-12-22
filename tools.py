import sys
import time


def humanize_time(tt):
    secs, microsecs = divmod(tt  ,  1)
    mins, secs      = divmod(secs, 60)
    hours, mins     = divmod(mins, 60)
    if hours:
        return '%02d:%02d:%02d.%02d' % (hours, mins, secs,100.0*microsecs)
    elif mins:
        return '%02d:%02d.%02d' % (mins, secs,100.0*microsecs)
    elif secs:
        return '%02d.%02d' % (secs,100.0*microsecs)
    else:
        return '%02d.%02d' % (0,100.0*microsecs)



class ProgressBar(object):
    """ProgressBar class holds the options of the progress bar.
    The options are:
        start   State from which start the progress. For example, if start is 
                5 and the end is 10, the progress of this state is 50%
        end     State in which the progress has terminated.
        width   --
        fill    String to use for "filled" used to represent the progress
        blank   String to use for "filled" used to represent remaining space.
        format  Format
        incremental
    """
    def __init__(self, start=0, end=10, width=12, fill='=', blank='.', format='[%(fill)s>%(blank)s] %(progress)s%%  Time Left: %(time)s', incremental=True):
        super(ProgressBar, self).__init__()

        self.start = start
        self.end = end
        self.width = width
        self.fill = fill
        self.blank = blank
        self.format = format
        self.incremental = incremental
        self.step = 100 / float(width) #fix
        self.reset()
        self.time = None
        self.init = True
        self.fin  = False
        self.t0   = None

    def __add__(self, increment):
        increment = self._get_progress(increment)
        if 100 > self.progress + increment:
            self.progress += increment
        else:
            self.progress = 100
            self.fin      = True
        return self

    def __str__(self):
        progressed = int(self.progress / self.step) #fix
        fill = progressed * self.fill
        blank = (self.width - progressed) * self.blank
        return self.format % {'fill': fill, 'blank': blank, 'progress': int(self.progress), 'time':str(self.time)}

    __repr__ = __str__

    def _get_progress(self, increment):
        return float(increment * 100) / self.end

    def reset(self):
        """Resets the current progress to the start point"""
        self.progress = self._get_progress(self.start)
        return self


class AnimatedProgressBar(ProgressBar):
    """Extends ProgressBar to allow you to use it straighforward on a script.
    Accepts an extra keyword argument named `stdout` (by default use sys.stdout)
    and may be any file-object to which send the progress status.
    """
    def __init__(self, *args, **kwargs):
        super(AnimatedProgressBar, self).__init__(*args, **kwargs)
        self.stdout = kwargs.get('stdout', sys.stdout)

    def show_progress(self):
        if self.init:
            self.t0   = time.time()
            self.init = False
            self.time = '-'
        else:
            self.time = humanize_time(((time.time() - self.t0)/self.progress)*(100.0-self.progress))
        if hasattr(self.stdout, 'isatty') and self.stdout.isatty():
            self.stdout.write('\r')
        else:
            self.stdout.write('\n')
        self.stdout.write(str(self))
        if self.fin:
            self.stdout.write('\n Total Time: '+humanize_time(time.time() - self.t0)+'\n')
        self.stdout.flush()


if __name__ == '__main__':
    p = AnimatedProgressBar(end=100, width=80)

    while True:
        p + 5
        p.show_progress()
        time.sleep(0.1)
        if p.progress == 100:
            break
    print #new line
