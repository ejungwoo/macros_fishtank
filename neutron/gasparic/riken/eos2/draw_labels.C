{

  int color=0;

  TText *p = new TText(30.,88.,"p");
  p->SetTextColor(color);
  p->Draw("same");
    
  TText *d = new TText(55.,95.,"d");
  d->SetTextColor(color);
  d->Draw("same");
    
  TText *t = new TText(80.,103.,"t");
  t->SetTextColor(color);
  t->Draw("same");
    
  TLatex *He4 = new TLatex(150.,78.,"^{4}He");
  He4->SetTextColor(color);
  He4->Draw("same");
    
  TLatex *He3 = new TLatex(110.,71.,"^{3}He");
  He3->SetTextColor(color);
  He3->Draw("same");
    
}
