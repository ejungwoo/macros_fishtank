{

  int color=0;

  TText *n = new TText(50.,60.,"n");
  n->SetTextColor(color);
  n->SetTextSize(0.07);
  n->Draw("same");
    
  TLatex *g = new TLatex(10.,35.,"#bf{#gamma}");
  g->SetTextColor(color);
  g->SetTextSize(0.08);
  g->Draw("same");

}
