
import sunpy
import os
from sunpy.net import hek
from sunpy.time import parse_time
import csv


def get_goes_event_list(tstart, tend, filename='goes_event_list.csv', goes_class_filter=None):
    """
    Function to get the GOES flare event list and save it as a csv.

    Parameters
    ----------
    tstart : ~str, datetime.datetime, astropy.time.Time
        start time of search
    tend : ~str, datetime.datetime, astropy.time.Time
        end time of search
    filename : ~str, optional
        name of csv file to save results to
    goes_class_filter : ~str, optional
        minimum GOES class flare to search for, only returns flares
        with a goes class greater than given

    Returns
    -------
    csv file

    Examples
    --------
    >>> t_start = '2014-01-01 00:00'
    >>> t_end = '2014-01-28 23:59'
    >>> goes_class_filter = 'C1.0'
    # searches for flares greater than C1.0 between t_start and t_end
    >>> get_goes_event_list(t_start, t_end, goes_class_filter=goes_class_filter)
    """

    client = hek.HEKClient()
    event_type = 'FL'

    if goes_class_filter:
        result = client.search(hek.attrs.Time(tstart, tend), hek.attrs.EventType(event_type), hek.attrs.FL.GOESCls >= goes_class_filter, hek.attrs.OBS.Observatory == 'GOES')
    else:
        result = client.search(hek.attrs.Time(tstart, tend), hek.attrs.EventType(event_type), hek.attrs.OBS.Observatory == 'GOES')

    goes_event_list = []

    for r in result:
        goes_event = {}
        tmpstring = str.split(str(r['event_starttime']), 'T')
        date = tmpstring[0]
        goes_event['event_date'] = date
        goes_event['start_time'] = parse_time(r['event_starttime'])
        goes_event['peak_time'] = parse_time(r['event_peaktime'])
        goes_event['end_time'] = parse_time(r['event_endtime'])
        goes_event['goes_class'] = str(r['fl_goescls'])
        goes_event_list.append(goes_event)

    with open(filename, 'w') as csvfile:
        fieldnames = goes_event_list[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(goes_event_list)):
            writer.writerow(goes_event_list[i])
