void add_merger_event(MyIDType Merging_ID1 , MyIDType Merging_ID2)
{
  memset(&MergerEvents[N_mergers], 0, sizeof(struct merger_data));
  MergerEvents[N_mergers].ID1 = Merging_ID1; 
  MergerEvents[N_mergers].ID2 = Merging_ID2;
  N_mergers++;
}
