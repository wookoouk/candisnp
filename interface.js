//define the style for the tooltip popover thingy

$(function () {
  $(".pop").popover(
  {
    trigger: 'hover',
    html: 'true',
    delay: 250,
    content: '<div style="font-size:6px;">'
  });
});

$(function () {
    $('.tip').tooltip()
});
