# svg sadness

Sadly, svg images are imperfectly rendered by some browsers.

Consider the following:

<img src="img/sadness.svg" alt="sadness" title="sadness" />

Gradually resize your window by making it narrower and narrower, making it as narrow as you can.
Watch the position of the right text end relative to the end of the surrounding rectangle.

For Firefox, the text always stays in the same relative position.

For Chrome, the text gets larger and larger relative to the rectangle.

For Safari, the relationship bounces around.

One would like to think that relative positions of objects within svg would maintain
a fixed relationship.  But this is not always the case.

This experiment was carried out on a Mac running High Sierra and might not apply to other platforms.

It is also possible that the imperfect rendering has something to do with how GitHub renders
markdown.
