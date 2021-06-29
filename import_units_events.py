from rec import *
from rec_db import *
from metadata import *

'''Prereqs
MetaData class and saved sessions, trials, events, units
generates conditions table with extra columns:
...
'''


def main():

    db = DataBase()

    sessions = db.tables['sessions']
    units = db.tables['units']
    activity = db.tables['activity']
    events = db.tables['events']

    for session in sessions.itertuples():

        
    [[timestamp_interval_within_activity(event_entry.StageTimeIndex,
                                         event_entry.StageTimeIndex + event_entry.StageDuration,
                                         db.tables['activity'].loc[unit_entry.Index])
      for event_entry in events.iloc.itertuples()] for unit_entry in units.iloc.itertuples()]

    # SAVE
    md.db_base_saver([('conditions', conditions)])


main()
