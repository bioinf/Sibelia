hll* merge2(hll* list1, hll* list2) {
  hll* totallist = 0;
  hll* temp;
  while (list1 || list2) {
    if ((list1 && !list2) || (list1->seq1start > list2->seq1start)) {
      temp = list1->next;
      list1->next = totallist;
      totallist = list1;
      list1 = temp;
    }
    else {
      temp = list2->next;
      list2->next = totallist;
      totallist = list2;
      list2 = temp;
    }
  }
  return totallist;
}

hll* findmiddle(hll* mylist) {
  hll* other = mylist;
  while (other && other->next) {
    other = other->next->next;
    mylist = mylist->next;
  }
  return mylist;
}

hll* sortList(hll* mylist) {
  hll* premid; 
  hll* mid;
  if (!mylist || !mylist->next)
    return mylist;

  premid = findmiddle(mylist);
  mid = premid->next;
  premid->next = 0;
  mylist = sortList(mylist);
  mid = sortList(mylist);
  mylist = merge2(mylist,mid);
}

