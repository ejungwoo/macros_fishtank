void sorthits()
{
  TObjArray array;

  STHit *hit;

  hit = new STHit(); hit -> SetS(1); array.Add(hit);
  hit = new STHit(); hit -> SetS(2); array.Add(hit);
  hit = new STHit(); hit -> SetS(3); array.Add(hit);
  hit = new STHit(); hit -> SetS(4); array.Add(hit);
  hit = new STHit(); hit -> SetS(2.2); array.Add(hit);
  hit = new STHit(); hit -> SetS(1.1); array.Add(hit);

  array.Sort();

  for (auto i = 0; i < 6; i++)
    cout << ((STHit *) array.At(i)) -> GetS() << endl;
}
