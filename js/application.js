$(function() {

  $('button.bifurcational').on('click', function(e) {
    e.preventDefault();
    $form = $(this).closest('form');

    params = getParams($form);
    
    calculateBifurcational(
      params['N'],
      params['R'],
      params['epsilon'],
      params['mu'],
      params['gamma'],
      params['phi'],
      params['x0'],
      params['x0_dot'],
      params['t0'],
      params['tn'],
      params['h'],
      params['lambda']
    );

  });

  $('button.method').on('click', function(e) {
    e.preventDefault();
    $form = $(this).closest('form');

    params = getParams($form);

    calculate(
      params['N'],
      params['p'],
      params['R'],
      params['epsilon'],
      params['mu'],
      params['gamma'],
      params['phi'],
      params['x0'],
      params['x0_dot'],
      params['t0'],
      params['tn'],
      params['h'],
      params['lambda']
    );
  });

});